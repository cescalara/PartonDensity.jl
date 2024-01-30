# # Fit with valence parametrisation
#
# In this example we show how to bring the PDF parametrisation and
# forward model together with `BAT.jl` to perform a fit of simulated data.

using BAT, DensityInterface
using PartonDensity
using QCDNUM
using Plots, Random, Distributions, ValueShapes, ParallelProcessingTools
using StatsBase, LinearAlgebra
include("../data/ZEUS_I1787035/ZEUS_I1787035.jl")
gr(fmt=:png);
rng = MersenneTwister(42)

# ## Simulate some data

# We can start off by simulating some fake data for us to fit. This way,
# we know exactly what initial conditions we have specified and can check
# the validity of our inference, assuming the generative model is the one that is producing our data.
#
# This is a good first check to work with.
#
# ### Specify the input PDFs
#
# See the *Input PDF parametrisation and priors* example for more information on the
# definition of the input PDFs. Here, we use the valence parametrisation.

weights = [5.0, 5.0, 1.0, 1.0, 1.0, 0.5, 0.5]
λ_u = 0.64;
K_u = 3.38;
λ_d = 0.67;
K_d = 4.73;
θ = get_θ_val(rng, λ_u, K_u, λ_d, K_d, weights)
pdf_params = ValencePDFParams(λ_u=λ_u, K_u=K_u, λ_d=λ_d, K_d=K_d,
    λ_g1=0.50, λ_g2=-0.63, K_g=4.23, λ_q=-0.23, K_q=5.0, θ=θ);

plot_input_pdfs(pdf_params)

# ### Go from PDFs to counts in ZEUS detector bins
#
# Given the input PDFs, we can then evolve, calculate the cross sections, and fold through
# the ZEUS transfer matrix to get counts in bins. Here, we make use of some simple helper
# functions to do so. For more details, see the *Forward model* example.

# first specify QCDNUM inputs
qcdnum_grid = QCDNUM.GridParams(x_min=[1.0e-3], x_weights=[1], nx=100,
    qq_bounds=[1.0e2, 3.0e4], qq_weights=[1.0, 1.0], nq=50, spline_interp=3)
qcdnum_params = QCDNUM.EvolutionParams(order=2, α_S=0.118, q0=100.0, grid_params=qcdnum_grid,
    n_fixed_flav=5, iqc=1, iqb=1, iqt=1, weight_type=1);

# now SPLINT and quark coefficients
splint_params = QCDNUM.SPLINTParams();
quark_coeffs = QuarkCoefficients();

# initialise QCDNUM
forward_model_init(qcdnum_params, splint_params)

# run forward model 
counts_pred_ep, counts_pred_em = forward_model(pdf_params, qcdnum_params, splint_params, quark_coeffs, MD_ZEUS_I1787035);

#
# take a poisson sample
nbins = size(counts_pred_ep)[1]
counts_obs_ep = zeros(UInt64, nbins)
counts_obs_em = zeros(UInt64, nbins)

for i in 1:nbins
    counts_obs_ep[i] = rand(rng, Poisson(counts_pred_ep[i]))
    counts_obs_em[i] = rand(rng, Poisson(counts_pred_em[i]))
end

#

plot(1:nbins, counts_pred_ep, label="Expected counts (eP)", color="blue")
plot!(1:nbins, counts_pred_em, label="Expected counts (eM)", color="red")
scatter!(1:nbins, counts_obs_ep, label="Detected counts (eP)", color="blue")
scatter!(1:nbins, counts_obs_em, label="Detected counts (eM)", color="red")
plot!(xlabel="Bin number")

# store
sim_data = Dict{String,Any}()
sim_data["nbins"] = nbins;
sim_data["counts_obs_ep"] = counts_obs_ep;
sim_data["counts_obs_em"] = counts_obs_em;

# write to file
pd_write_sim("output/simulation.h5", pdf_params, sim_data)
QCDNUM.save_params("output/params_val.h5", qcdnum_params)
QCDNUM.save_params("output/params_val.h5", splint_params)

# ## Fit the simulated data
#
# Now we can try to fit this simulated data using `Bat.jl`.
# The first step is to define the prior and likelihood.
# For now, let's try relatively narrow priors centred on the true values.

prior = NamedTupleDist(
    θ_tmp=Dirichlet(weights),
    λ_u=Truncated(Normal(pdf_params.λ_u, 1), 0, 1),
    K_u=Truncated(Normal(pdf_params.K_u, 1), 2, 10),
    λ_d=Truncated(Normal(pdf_params.λ_d, 1), 0, 1),
    K_d=Truncated(Normal(pdf_params.K_d, 1), 2, 10),
    λ_g1=Truncated(Normal(pdf_params.λ_g1, 1), -1, 0),
    λ_g2=Truncated(Normal(pdf_params.λ_g2, 1), -1, 0),
    K_g=Truncated(Normal(pdf_params.K_g, 1), 2, 10),
    λ_q=Truncated(Normal(pdf_params.λ_q, 0.1), -1, 0),
    K_q=Truncated(Normal(pdf_params.K_q, 0.5), 3, 7)
);

# The likelihood is similar to that used in the *input PDF parametrisation* example.
# We start by accessing the current parameter set of the sampler's iteration,
# then running the forward model to get the predicted counts and comparing to
# the observed counts using a simple Poisson likelihood.
#

likelihood = let d = sim_data

    counts_obs_ep = d["counts_obs_ep"]
    counts_obs_em = d["counts_obs_em"]
    nbins = d["nbins"]

    logfuncdensity(function (params)

        θ = get_scaled_θ(params.λ_u, params.K_u, params.λ_d,
            params.K_d, Vector(params.θ_tmp))

        pdf_params = ValencePDFParams(λ_u=params.λ_u, K_u=params.K_u, λ_d=params.λ_d,
            K_d=params.K_d, λ_g1=params.λ_g1, λ_g2=params.λ_g2,
            K_g=params.K_g, λ_q=params.λ_q, K_q=params.K_q, θ=θ)

        counts_pred_ep, counts_pred_em = @critical forward_model(pdf_params, qcdnum_params, splint_params, quark_coeffs, MD_ZEUS_I1787035)

        ll_value = 0.0
        for i in 1:nbins
            ll_value += logpdf(Poisson(counts_pred_ep[i]), counts_obs_ep[i])
            ll_value += logpdf(Poisson(counts_pred_em[i]), counts_obs_em[i])
        end

        return ll_value
    end)
end

# We can now run the MCMC sampler. We will start by using the
# Metropolis-Hastings algorithm as implemented in `BAT.jl`.
# To get reasonable results, we need to run the sampler for a
# long time (several hours). To save time in this demo, we will
# work with a ready-made results file. To actually run the sampler,
# simply uncomment the code below.

#posterior = PosteriorDensity(likelihood, prior);

#samples = bat_sample(posterior, MCMCSampling(mcalg=MetropolisHastings(), nsteps=10^5, nchains=2)).result;

# Alternatively, we could also try a nested sampling approach
# here for comparison. This is easily done thanks to the
# interface of `BAT.jl`, you will just need to add the
# `NestedSamplers.jl` package.

#import NestedSamplers
#samples = bat_sample(posterior, EllipsoidalNestedSampling()).result

# If you run the sampler, be sure to save
# the results for further analysis

#import HDF5
#bat_write("output/results.h5", samples)

# ## Analysing the results
#
# First, let's load our simulation inputs and results

pdf_params, sim_data = pd_read_sim("output/demo_simulation_valence.h5");
samples = bat_read("output/demo_results_valence.h5").result;

# We can use the same QCDNUM params as above
loaded_params = QCDNUM.load_params("output/params_val.h5")
qcdnum_params = loaded_params["evolution_params"]
splint_params = loaded_params["splint_params"]

# We can check some diagnostics using built in `BAT.jl`, such as the
# effective sample size shown below

bat_eff_sample_size(unshaped.(samples))[1]

# We see a value for each of our 15 total parameters. As the
# Metropolis-Hastings algorithm's default implementation
# isn't very efficient, we see that the effective sample size
# is only a small percentage of the input `nsteps`. We should try
# to improve this if possible, or use a much larger `nsteps` value.
#
# For demonstration purposes, we will continue to show how we can
# visualise the results in this case. For robust inference, we need
# to improve the sampling stage above.

# We can use `BAT.jl`'s built in plotting recipes to show the marginals,
# for example, consider `λ_u`, and compare to the known truth.

plot(
    samples, :(λ_u),
    nbins=50,
    colors=[:skyblue4, :skyblue3, :skyblue1],
    alpha=0.7,
    marginalmode=false,
    legend=:topleft
)
vline!([pdf_params.λ_u], color="black", label="truth", lw=3)

# If we want to compare the momentum weights, we must transform
# from `θ_tmp` to `θ`, as shown below. Here, we transform using
# a helper function, convert the result to a matrix, and
# access the ith weight with the integer `i`.

θ = [get_scaled_θ(λ_u, K_u, λ_d, K_d, Vector(θ_tmp)) for (λ_u, K_u, λ_d, K_d, θ_tmp)
     in
     zip(samples.v.λ_u, samples.v.K_u, samples.v.λ_d, samples.v.K_d, samples.v.θ_tmp)];

θ = transpose(reduce(vcat, transpose.(θ)))

i = 1
hist = append!(Histogram(0:0.02:1), θ[i, :])
plot(
    normalize(hist, mode=:density),
    st=:steps, label="Marginal posterior"
)
vline!([pdf_params.θ[i]], color="black", label="truth", lw=3)

# Rather than making a large plot 15 different marginals,
# it can be more useful to visualise the posterior distribution
# in differently, such as the shape of the distributions
# we are trying to fit, or the *model space*. Helper functions
# exist for doing just this.

# Using BAT recipe
function wrap_xtotx(p::NamedTuple, x::Real)
    θ = get_scaled_θ(p.λ_u, p.K_u, p.λ_d, p.K_d, Vector(p.θ_tmp))
    pdf_params = ValencePDFParams(λ_u=p.λ_u, K_u=p.K_u, λ_d=p.λ_d, K_d=p.K_d, λ_g1=p.λ_g1,
        λ_g2=p.λ_g2, K_g=p.K_g, λ_q=p.λ_q, K_q=p.K_q, θ=θ)
    return log(xtotx(x, pdf_params))
end

x_grid = range(1e-3, stop=1, length=50)
plot(x_grid, wrap_xtotx, samples, colors=[:skyblue4, :skyblue3, :skyblue1],
    legend=:topright)
plot!(x_grid, [log(xtotx(x, pdf_params)) for x in x_grid], color="black", lw=3,
    label="Truth", linestyle=:dash)
plot!(ylabel="log(xtotx)")

# Using `PartonDensity.jl`
plot_model_space(pdf_params, samples, nsamples=500)

# Alternatively, we can also visualise the implications of the fit
# in the *data space*, as shown below. 

plot_data_space(pdf_params, sim_data, samples, qcdnum_params,
    splint_params, quark_coeffs, nsamples=500)
