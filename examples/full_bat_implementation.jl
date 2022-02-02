# # A fit with BAT.jl
#
# In this example we show how to bring the PDF parametrisation and
# forward model together with `BAT.jl` to perform a fit of simulated data.
# This fit is a work in progress and just a starting point for verification
# of the method.

using BAT, DensityInterface
using PartonDensity
using QCDNUM
using Plots, Random, Distributions, ValueShapes, ParallelProcessingTools

# ## Simulate some data

seed = 42
Random.seed!(seed) # for reproducibility

# We can start off by simulating some fake data for us to fit. This way,
# we know exactly what initial conditions we have specified and can check
# the validity of our inference, assuming the generative model is the one that is producing our data.
#
# This is a good first check to work with.
#
# ### Specify the input PDFs
#
# See the *Input PDF parametrisation and priors* example for more information on the
# definition of the input PDFs.

pdf_params = PDFParameters(λ_u=0.7, K_u=4.0, λ_d=0.5, K_d=6.0,
    λ_g1=0.7, λ_g2=-0.4, K_g=6.0, λ_q=-0.5, weights=[1, 0.5, 0.3, 0.2, 0.1, 0.1, 0.1]);

#plot_input_pdfs(pdf_params)

# ### Go from PDFs to counts in ZEUS detector bins
#
# Given the input PDFs, we can then evolve, calculate the cross sections, and fold through
# the ZEUS transfer matrix to get counts in bins. Here, we make use of some simple helper
# functions to do so. For more details, see the *Forward model* example.

# first specify QCDNUM inputs
qcdnum_grid = QCDNUMGrid(x_min=[1.0e-3], x_weights=[1], nx=100,
    qq_bounds=[1.0e2, 3.0e4], qq_weights=[1.0, 1.0], nq=50, spline_interp=3)
qcdnum_params = QCDNUMParameters(order=2, α_S=0.118, q0=100.0, grid=qcdnum_grid,
    n_fixed_flav=5, iqc=1, iqb=1, iqt=1, weight_type=1);

# now SPLINT and quark coefficients
splint_params = SPLINTParameters();
quark_coeffs = QuarkCoefficients();

# initialise QCDNUM
forward_model_init(qcdnum_grid, qcdnum_params, splint_params)

# run forward model 
counts_pred_ep, counts_pred_em = forward_model(pdf_params, qcdnum_params, 
    splint_params, quark_coeffs);

#
# take a poisson sample
nbins = size(counts_pred_ep)[1]
counts_obs_ep = zeros(UInt64, nbins)
counts_obs_em = zeros(UInt64, nbins)

for i in 1:nbins
    counts_obs_ep[i] = rand(Poisson(counts_pred_ep[i]))
    counts_obs_em[i] = rand(Poisson(counts_pred_em[i]))
end

#

plot(1:nbins, counts_pred_ep, label="Expected counts (eP)", color="blue")
plot!(1:nbins, counts_pred_em, label="Expected counts (eM)", color="red")
scatter!(1:nbins, counts_obs_ep, label="Detected counts (eP)", color="blue")
scatter!(1:nbins, counts_obs_em, label="Detected counts (eM)", color="red")
plot!(xlabel="Bin number")

# store
sim_data = Dict{String, Any}()
sim_data["nbins"] = nbins;
sim_data["counts_obs_ep"] = counts_obs_ep;
sim_data["counts_obs_em"] = counts_obs_em;

# write to file
pd_write_sim("output/simulation.h5", pdf_params, sim_data)

# ## Fit the simulated data
#
# Now we can try to fit this simulated data using `Bat.jl`.
# The first step is to define the prior and likelihood.
# For now, let's try relatively narrow priors centred on the true values.

prior = NamedTupleDist(
    θ = Dirichlet([1, 0.5, 0.3, 0.2, 0.1, 0.1, 0.1]),
    λ_u = Truncated(Normal(pdf_params.λ_u, 0.2), 0, 1), 
    K_u = Truncated(Normal(pdf_params.K_u, 1), 2, 10),
    λ_d = Truncated(Normal(pdf_params.λ_d, 0.2), 0, 1), 
    K_d = Truncated(Normal(pdf_params.K_d, 1), 2, 10),
    λ_g1 = Truncated(Normal(pdf_params.λ_g1, 0.2), 0, 1),
    λ_g2 = Truncated(Normal(pdf_params.λ_g2, 0.2), -1, 0),
    K_g =  Truncated(Normal(pdf_params.K_g, 1), 2, 10),
    λ_q = Truncated(Normal(pdf_params.λ_q, 0.1), -1, 0),
);

# The likelihood is similar to that used in the *input PDF parametrisation* example.
# We start by accessing the current parameter set of the sampler's iteration,
# then running the forward model to get the predicted counts and comparing to
# the observed counts using a simple Poisson likelihood.
#
# The `@critical` macro is used because `forward_model()` is currently not thread safe, so
# this protects it from being run in parallel.

likelihood = let d = sim_data

    counts_obs_ep = d["counts_obs_ep"]
    counts_obs_em = d["counts_obs_em"]
    nbins = d["nbins"]

    logfuncdensity(function (params)

            pdf_params = PDFParameters(λ_u=params.λ_u, K_u=params.K_u, λ_d=params.λ_d,
                                       K_d=params.K_d, λ_g1=params.λ_g1, λ_g2=params.λ_g2,
                                       K_g=params.K_g, λ_q=params.λ_q, θ=Vector(params.θ))

            @critical counts_pred_ep, counts_pred_em = forward_model(pdf_params, 
                qcdnum_params, splint_params, quark_coeffs);

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
# To demonstrate, we run a short chain of 10^3 iterations, which
# should take around 10 minutes.

posterior = PosteriorDensity(likelihood, prior);
samples = bat_sample(
    posterior,
    MCMCSampling(mcalg=MetropolisHastings(), nsteps=10^3, nchains=2)
).result;

# Let's save the result for further analysis

import HDF5
bat_write("output/results.h5", samples)

# ## Analysing the results
#
# First, let's load our simulation inputs and results

pdf_params, sim_data = pd_read_sim("output/simulation.h5");
samples = bat_read("output/results.h5").result;

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

# Rather than making a large plot 15 different marginals,
# it can be more useful to visualise the posterior distribution
# in differently, such as the shape of the distributions
# we are trying to fit, or the *model space*. Helper functions
# exist for doing just this.

# Using BAT recipe
function wrap_xtotx(p::NamedTuple{(:K_d, :K_g, :K_u, :θ, :λ_d, :λ_g1, 
                                   :λ_g2, :λ_q, :λ_u)}, x::Real)
    pdf_params = PDFParameters(λ_u=p.λ_u, K_u=p.K_u, λ_d=p.λ_d, K_d=p.K_d, λ_g1=p.λ_g1, 
        λ_g2=p.λ_g2, K_g=p.K_g, λ_q=p.λ_q, θ=p.θ)
    return xtotx(x, pdf_params)
end

x_grid = range(1e-3, stop=1, length=50)
plot(x_grid, wrap_xtotx, samples, colors=[:skyblue4, :skyblue3, :skyblue1], 
     legend=:topright)
plot!(x_grid, [xtotx(x, pdf_params) for x in x_grid], color="black", lw=3, 
      label="Truth", linestyle=:dash)
plot!(xaxis=:log)

# Using `PartonDensity.jl`
plot_model_space(pdf_params, samples)

# Alternatively, we can also visualise the implications of the fit
# in the *data space*, as shown below. 

plot_data_space(pdf_params, sim_data, samples, qcdnum_grid, qcdnum_params, 
                splint_params, quark_coeffs)

# The first results seem promising, but these are really just first checks
# and more work will have to be done to verify the method.
