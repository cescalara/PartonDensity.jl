using BAT, DensityInterface
using PartonDensity
using QCDNUM
using Plots, Random, Distributions, ValueShapes, ParallelProcessingTools
using StatsBase, LinearAlgebra

gr(fmt=:png);

seed = 42
Random.seed!(seed) # for reproducibility

#pdf_params = ValencePDFParams(λ_u=0.64, K_u=3.38, λ_d=0.67, K_d=4.73,
#                              λ_g1=0.50, λ_g2=-0.63, K_g=4.23, λ_q=-0.23, weights=[5., 5., 1., 1., 1., 0.5, 0.5]);

#plot_input_pdfs(pdf_params)
pdf_params = DirichletPDFParams(K_u=4.0, K_d=4.0, λ_g1=1.5, λ_g2=-0.4, K_g=6.0, 
                                λ_q=-0.25, weights=[30., 15., 12., 6., 3.6, 0.85, 0.85, 0.85, 0.85]);

@info "Valence λ:" pdf_params.λ_u pdf_params.λ_d


#qcdnum_grid = QCDNUMGrid(x_min=[1.0e-3], x_weights=[1], nx=100,
#    qq_bounds=[1.0e2, 3.0e4], qq_weights=[1.0, 1.0], nq=50, spline_interp=3)
#qcdnum_params = QCDNUMParameters(order=2, α_S=0.118, q0=100.0, grid=qcdnum_grid,
#    n_fixed_flav=5, iqc=1, iqb=1, iqt=1, weight_type=1);
# first specify QCDNUM inputs
qcdnum_grid = QCDNUMGrid(x_min=[1.0e-3, 1.0e-1, 5.0e-1], x_weights=[1, 2, 2], nx=100,
                         qq_bounds=[1.0e2, 3.0e4], qq_weights=[1.0, 1.0], nq=50, spline_interp=3)
qcdnum_params = QCDNUMParameters(order=2, α_S=0.118, q0=100.0, grid=qcdnum_grid,
                                 n_fixed_flav=5, iqc=1, iqb=1, iqt=1, weight_type=1);


splint_params = SPLINTParameters();
quark_coeffs = QuarkCoefficients();

forward_model_init_sysErr(qcdnum_grid, qcdnum_params, splint_params)

counts_pred_ep, counts_pred_em = forward_model(pdf_params, qcdnum_params,
    splint_params, quark_coeffs);

nbins = size(counts_pred_ep)[1]
counts_obs_ep = zeros(UInt64, nbins)
counts_obs_em = zeros(UInt64, nbins)

for i in 1:nbins
    counts_obs_ep[i] = rand(Poisson(counts_pred_ep[i]))
    counts_obs_em[i] = rand(Poisson(counts_pred_em[i]))
end

plot(1:nbins, counts_pred_ep, label="Expected counts (eP)", color="blue")
plot!(1:nbins, counts_pred_em, label="Expected counts (eM)", color="red")
scatter!(1:nbins, counts_obs_ep, label="Detected counts (eP)", color="blue")
scatter!(1:nbins, counts_obs_em, label="Detected counts (eM)", color="red")
plot!(xlabel="Bin number")

sim_data = Dict{String, Any}()
sim_data["nbins"] = nbins;
sim_data["counts_obs_ep"] = counts_obs_ep;
sim_data["counts_obs_em"] = counts_obs_em;

#pd_write_sim("simulation.h5", pdf_params, sim_data)

prior = NamedTupleDist(
    θ_tmp = Dirichlet(pdf_params.weights),
    λ_u = Truncated(Normal(pdf_params.λ_u, 1), 0, 1),
    K_u = Truncated(Normal(pdf_params.K_u, 1), 2, 10),
    λ_d = Truncated(Normal(pdf_params.λ_d, 1), 0, 1),
    K_d = Truncated(Normal(pdf_params.K_d, 1), 2, 10),
    λ_g1 = Truncated(Normal(pdf_params.λ_g1, 1), -1, 0),
    λ_g2 = Truncated(Normal(pdf_params.λ_g2, 1), -1, 0),
    K_g =  Truncated(Normal(pdf_params.K_g, 1), 2, 10),
    λ_q = Truncated(Normal(pdf_params.λ_q, 0.1), -1, 0),
    beta0_1=  Truncated(Normal(0, 1), -5, 5), 
    beta0_2=   Truncated(Normal(0, 1), -5, 5),    
    beta0_3= Truncated(Normal(0, 1), -5, 5), 
    beta0_4=  Truncated(Normal(0, 1), -5, 5), 
    beta0_5=  Truncated(Normal(0, 1), -5, 5), 
    beta0_6=  Truncated(Normal(0, 1), -5, 5), 
    beta0_7=  Truncated(Normal(0, 1), -5, 5), 
    beta0_8=   Truncated(Normal(0, 1), -5, 5)
);

likelihood = let d = sim_data

    counts_obs_ep = d["counts_obs_ep"]
    counts_obs_em = d["counts_obs_em"]
    nbins = d["nbins"]

    logfuncdensity(function (params)

            θ = get_scaled_θ(params.λ_u, params.K_u, params.λ_d,
                             params.K_d, Vector(params.θ_tmp))

            pdf_params = ValencePDFParams(λ_u=params.λ_u, K_u=params.K_u, λ_d=params.λ_d,
                                          K_d=params.K_d, λ_g1=params.λ_g1, λ_g2=params.λ_g2,
                                          K_g=params.K_g, λ_q=params.λ_q, θ=θ)

            
           # ParErrs = zeros(8)
        ParErrs = [params.beta0_1,params.beta0_2,params.beta0_3,params.beta0_4,
                params.beta0_5,params.beta0_6,params.beta0_7,params.beta0_8]
            
            counts_pred_ep, counts_pred_em = @critical forward_model_sysErr(pdf_params,
                qcdnum_params, splint_params, quark_coeffs,ParErrs );

            ll_value = 0.0
            for i in 1:nbins

                if counts_pred_ep[i] < 0
                   #@warn "counts_pred_ep[i] < 0, setting to 0" i counts_pred_ep[i]
                   counts_pred_ep[i] = 0
                end

                if counts_pred_em[i] < 0
                   #@warn "counts_pred_em[i] < 0, setting to 0" i counts_pred_em[i]
                   counts_pred_em[i] = 0
                end

                ll_value += logpdf(Poisson(counts_pred_ep[i]), counts_obs_ep[i])
                ll_value += logpdf(Poisson(counts_pred_em[i]), counts_obs_em[i])
            end

            return ll_value
    end)
end

posterior = PosteriorDensity(likelihood, prior);
convergence = BrooksGelmanConvergence(threshold=1.3);
burnin = MCMCMultiCycleBurnin(max_ncycles=100);

samples = bat_sample(posterior, MCMCSampling(mcalg=MetropolisHastings(), nsteps=10^4, nchains=2)).result;

#import NestedSamplers
#samples = bat_sample(posterior, EllipsoidalNestedSampling()).result

#import HDF5
#bat_write("output/results.h5", samples)

#pdf_params, sim_data = pd_read_sim("output/demo_simulation.h5");
#samples = bat_read("output/demo_results.h5").result;

#bat_eff_sample_size(unshaped.(samples))[1]

# plot(
#     samples, :(λ_u),
#     nbins=50,
#     colors=[:skyblue4, :skyblue3, :skyblue1],
#     alpha=0.7,
# #     marginalmode=false,
#     legend=:topleft
# )
# vline!([pdf_params.λ_u], color="black", label="truth", lw=3)

# θ = [get_scaled_θ(λ_u, K_u, λ_d, K_d, Vector(θ_tmp)) for (λ_u, K_u, λ_d, K_d, θ_tmp)
#      in zip(samples.v.λ_u, samples.v.K_u, samples.v.λ_d, samples.v.K_d, samples.v.θ_tmp)];

# θ = transpose(reduce(vcat,transpose.(θ)))

# i = 1
# hist = append!(Histogram(0:0.02:1), θ[i,:])
# plot(
#     normalize(hist, mode=:density),
#     st = :steps, label="Marginal posterior"
# )
# vline!([pdf_params.θ[i]], color="black", label="truth", lw=3)

# function wrap_xtotx(p::NamedTuple, x::Real)
#     θ = get_scaled_θ(p.λ_u, p.K_u, p.λ_d, p.K_d, Vector(p.θ_tmp))
#     pdf_params = ValencePDFParams(λ_u=p.λ_u, K_u=p.K_u, λ_d=p.λ_d, K_d=p.K_d, λ_g1=p.λ_g1,
#         λ_g2=p.λ_g2, K_g=p.K_g, λ_q=p.λ_q, θ=θ)
#     return log(xtotx(x, pdf_params))
# end

# x_grid = range(1e-3, stop=1, length=50)
# plot(x_grid, wrap_xtotx, samples, colors=[:skyblue4, :skyblue3, :skyblue1],
#      legend=:topright)
# plot!(x_grid, [log(xtotx(x, pdf_params)) for x in x_grid], color="black", lw=3,
#       label="Truth", linestyle=:dash)
# plot!(ylabel="log(xtotx)")

#plot_model_space(pdf_params, samples, nsamples=500)

#plot_data_space(pdf_params, sim_data, samples, qcdnum_grid, qcdnum_params,
#                splint_params, quark_coeffs, nsamples=500)