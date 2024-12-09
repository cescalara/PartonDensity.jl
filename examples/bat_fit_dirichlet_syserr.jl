# # Fit with Dirichlet parametrisation and systematic errors
#
# In this example we show how to bring the PDF parametrisation and
# forward model together with `BAT.jl` to perform a fit of simulated data.
#
# This time, we also include the effect of systematic uncertainties. 
# This is done as described in [arXiv:2209.06571](https://arxiv.org/abs/2209.06571).

using PartonDensity
using BAT, DensityInterface
using QCDNUM
using Plots, Random, Distributions, ValueShapes, ParallelProcessingTools
using StatsBase, LinearAlgebra

zeus_include_path = string(chop(pathof(PartonDensity), tail=20), "data/ZEUS_I1787035/ZEUS_I1787035.jl")
MD_ZEUS_I1787035 = include(zeus_include_path)
gr(fmt=:png);
rng = MersenneTwister(42);

weights = [30.0, 15.0, 12.0, 6.0, 3.6, 0.85, 0.85, 0.85, 0.85]
θ = rand(rng, Dirichlet(weights))
pdf_params = DirichletPDFParams(K_u=4.0, K_d=4.0, λ_g1=1.5, λ_g2=-0.4, K_g=6.0,
    λ_q=-0.25, K_q=5.0, θ=θ);

@info "Valence λ:" pdf_params.λ_u pdf_params.λ_d

# first specify QCDNUM inputs
qcdnum_grid = QCDNUM.GridParams(x_min=[1.0e-3, 1.0e-1, 5.0e-1], x_weights=[1, 2, 2], nx=100,
    qq_bounds=[1.0e2, 3.0e4], qq_weights=[1.0, 1.0], nq=50, spline_interp=3)
qcdnum_params = QCDNUM.EvolutionParams(order=2, α_S=0.118, q0=100.0, grid_params=qcdnum_grid,
    n_fixed_flav=5, iqc=1, iqb=1, iqt=1, weight_type=1);

splint_params = QCDNUM.SPLINTParams();
quark_coeffs = QuarkCoefficients();

# We include the effects of systematic errors into the simulation, by sampling 
# from a multivariate normal with mean 0 and standard deviation 1.
#
# These factors are then applied to a precomputed matrix in order to scale the 
# expected counts accordingly.
forward_model_init(qcdnum_params, splint_params)

sys_err_params = rand(rng, MvNormal(zeros(nsyst), ones(nsyst)))
counts_pred_ep, counts_pred_em = forward_model(pdf_params, qcdnum_params, splint_params, quark_coeffs, MD_ZEUS_I1787035, sys_err_params);

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

sim_data = Dict{String,Any}()
sim_data["nbins"] = nbins;
sim_data["counts_obs_ep"] = counts_obs_ep;
sim_data["counts_obs_em"] = counts_obs_em;

pd_write_sim("output/simulation.h5", pdf_params, sim_data, MD_ZEUS_I1787035)
QCDNUM.save_params("output/params_dir_sys.h5", qcdnum_params)
QCDNUM.save_params("output/params_dir_sys.h5", splint_params)

# Here, we include a prior over the 8 systematic error parameters, such that we marginalise over them.
prior = NamedTupleDist(
    θ=Dirichlet(weights),
    K_u=Uniform(3.0, 7.0),
    K_d=Uniform(3.0, 7.0),
    λ_g1=Uniform(1.0, 2.0),
    λ_g2=Uniform(-0.5, -0.1),
    K_g=Uniform(3.0, 7.0),
    λ_q=Uniform(-0.5, -0.1),
    K_q=Uniform(3.0, 7.0),
    beta0_1=Truncated(Normal(0, 1), -5, 5),
    beta0_2=Truncated(Normal(0, 1), -5, 5),
    beta0_3=Truncated(Normal(0, 1), -5, 5),
    beta0_4=Truncated(Normal(0, 1), -5, 5),
    beta0_5=Truncated(Normal(0, 1), -5, 5),
    beta0_6=Truncated(Normal(0, 1), -5, 5),
    beta0_7=Truncated(Normal(0, 1), -5, 5),
    beta0_8=Truncated(Normal(0, 1), -5, 5)
);

likelihood = let d = sim_data

    counts_obs_ep = d["counts_obs_ep"]
    counts_obs_em = d["counts_obs_em"]
    nbins = d["nbins"]

    logfuncdensity(function (params)

        pdf_params = DirichletPDFParams(K_u=params.K_u, K_d=params.K_d, λ_g1=params.λ_g1, λ_g2=params.λ_g2,
            K_g=params.K_g, λ_q=params.λ_q, K_q=params.K_q, θ=Vector(params.θ))

        #The sys_err_params must also be passed to the forward model here.
        sys_err_params = [params.beta0_1, params.beta0_2, params.beta0_3, params.beta0_4,
            params.beta0_5, params.beta0_6, params.beta0_7, params.beta0_8]

        counts_pred_ep, counts_pred_em = @critical forward_model(pdf_params, qcdnum_params, splint_params, quark_coeffs, MD_ZEUS_I1787035, sys_err_params)

        ll_value = 0.0
        for i in 1:nbins

            if counts_pred_ep[i] < 0
                @debug "counts_pred_ep[i] < 0, setting to 0" i counts_pred_ep[i]
                counts_pred_ep[i] = 0
            end
            if counts_pred_em[i] < 0
                @debug "counts_pred_em[i] < 0, setting to 0" i counts_pred_em[i]
                counts_pred_em[i] = 0
            end

            ll_value += logpdf(Poisson(counts_pred_ep[i]), counts_obs_ep[i])
            ll_value += logpdf(Poisson(counts_pred_em[i]), counts_obs_em[i])
        end

        return ll_value
    end)
end

# Check that we can evaluate the posterior:
posterior = lbqintegral(likelihood, prior)
BAT.checked_logdensityof(posterior, rand(prior))

# We can now run the MCMC sampler. We will start by using the
# Metropolis-Hastings algorithm as implemented in `BAT.jl`.
# To get reasonable results, we need to run the sampler for a
# long time (several hours). To actually run the sampler,
# simply uncomment the code below. To see how to work with 
# demo output results, check out the other fit examples.

#mcalg = MetropolisHastings()
#convergence = BrooksGelmanConvergence(threshold=1.3)
#samples = bat_sample(posterior, MCMCSampling(mcalg=mcalg, nsteps=100, nchains=2, strict=false)).result;

#import HDF5
#bat_write("output/results.h5", samples)

# Examples of how to analyse the results of the fit can be found in the 
# *Fit with valence parametrisation* and  *Fit with Dirichlet parametrisation* 
# sections of this documentation