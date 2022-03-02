# # Prior model implementation in BAT.jl
#
# In this notebook, the implementation of the prior model is demonstrated,
# starting with a simple two-component model and scaling up to the full 9 component model.

using Distributions, StatsBase, LinearAlgebra
using Plots, SpecialFunctions, Printf, Random, ValueShapes
using BAT, DensityInterface
const sf = SpecialFunctions;

gr(fmt=:png); 

# ## Simple two-component model: gluons 
#
# First we start with a simpler two-component model, and show all the
# steps explicity for clarity.
#
# ### Forward model
#
# The gluon distributions are parametrised by
# ```math
# x g(x) = A_{g1} x^{\lambda_g1}(1-x)^{K_g} + A_{g2} x^{\lambda_{g2}}(1-x)^5.
# ```
#
# We also want to impose 
#
# ```math
# \int_0^1 x g(x) dx = A_{g1} B(\lambda_{g1}+1, K_g+1) + A_{g2} B(\lambda_{g2}+1, 1+5) = 1,
# ```
#
# with B(.,.) the Beta function.
#
# Start by defining some useful functions:

function xg1x(x, λ_g1, K_g, θ_1)
    A_g1 = θ_1 / sf.beta(λ_g1 + 1, K_g + 1)
    return A_g1 * x^λ_g1 * (1 - x)^K_g
end

function xg2x(x, λ_g2, θ_2)
    A_g2 = θ_2 / sf.beta(λ_g2 + 1, 5 + 1)
    return A_g2 * x^λ_g2 * (1 - x)^5
end

function xgx(x, λ_g1, λ_g2, K_g, θ)
    xg1 = xg1x(x, λ_g1, K_g, θ[1])
    xg2 = xg2x(x, λ_g2, θ[2])
    return xg1 + xg2
end

# Choose true values for the high-level parameters and show what the resulting model looks like.

θ = [0.5, 0.5]
λ_g1 = 0.5 # rand(Uniform(0, 1))
λ_g2 = -0.7 # rand(Uniform(-1, 0))
K_g = 3 # rand(Uniform(2, 10))
truths = (θ = θ, λ_g1 = λ_g1, λ_g2 = λ_g2, K_g = K_g);

A_g1 = θ[1] / sf.beta(λ_g1+1, K_g+1)
A_g2 = θ[2] / sf.beta(λ_g2+1, 5+1);

# Check integral = 1
total = A_g1 * sf.beta(λ_g1+1, K_g+1) + A_g2 * sf.beta(λ_g2+1, 5+1)
print("Integral = ", total)

# Plot true model
x_grid = range(0, stop=1, length=50)

xg1 = A_g1 * x_grid.^λ_g1 .* (1 .- x_grid).^K_g
xg2 = A_g2 * x_grid.^λ_g2 .* (1 .- x_grid).^5

plot(x_grid, [xg1x(x, λ_g1, K_g, θ[1]) for x in x_grid], 
      alpha=0.7, label="x g1(x)", lw=3, color="green")
plot!(x_grid, [xg2x(x, λ_g2, θ[2]) for x in x_grid], 
     alpha=0.7, label="x g2(x)", lw=3, color="blue")
plot!(x_grid, [xgx(x, λ_g1, λ_g2, K_g, θ) for x in x_grid],
      alpha=0.7, label="x g1(x) + x g2(x)", lw=3, color="red")
plot!(xlabel="x")

# Now, for the purposes of testing the prior implementation,
# sample some data from this distribution assuming that the data
# are produced by integrating over the function in different bins,
# and multiplying by some factor.
# Then, plot the model and data to compare.

seed = 42
Random.seed!(seed) # for reproducibility

bins = 0.0:0.05:1.0
bin_widths = bins[2:end] - bins[1:end-1]
bin_centers = (bins[1:end-1] + bins[2:end]) / 2

N = 1000
nbins = size(bin_centers)[1]

expected_counts = zeros(nbins)
observed_counts = zeros(Integer, nbins)
for i in 1:nbins
    xg = xgx(bin_centers[i], λ_g1, λ_g2, K_g, θ) * N
    expected_counts[i] = bin_widths[i] * xg 
    observed_counts[i] = rand(Poisson(expected_counts[i]))
end

plot(bin_centers, [xgx(x, λ_g1, λ_g2, K_g, θ) for x in bin_centers] .* bin_widths * N,
    alpha=0.7, label="Expected", lw=3, color="red")
scatter!(bin_centers, observed_counts, lw=3, label="Observed", color="black")

# Store the data in a simple dict to pass to the likelihood later.

data = Dict()
data["N"] = N
data["bin_centers"] = bin_centers;
data["observed_counts"] = observed_counts;
data["bin_widths"] = bin_widths;

# ### Fit
#
# To fit this example data, we choose a prior over our hyperparameters `θ`, `λ_g1`, `λ_g2` and `K_g`.
#
# We decide to choose a sensible Dirichlet prior, and have a look at
# some samples to help understand what this means.

dirichlet = Dirichlet([1, 1])
test = rand(dirichlet, 1000)
plot(append!(Histogram(0:0.1:1), test[1, :]))
plot!(append!(Histogram(0:0.1:1), test[2, :]))

# So, we see that this means all weights are equally likely for
# both components.
#
# Now we can define the prior:
#
# ```math
# \mathrm{Prior} - P(\theta) P(\lambda_{g1}) P(\lambda_{g2}) P(K_g).
# ```

prior = NamedTupleDist(
    θ = Dirichlet([1, 1]),
    λ_g1 = Uniform(0, 1), #Truncated(Normal(0.7, 0.1), 0, 1)
    λ_g2 = Uniform(-1, 0), #Truncated(Normal(-0.7, 0.1), -1, 0),
    K_g =  Uniform(2, 10), # Truncated(Normal(1, 2), 2, 10),
);

# Define a simple poisson likelihood that describes the test data that we generated above.

likelihood = let d = data, f = xgx

    observed_counts = d["observed_counts"]
    bin_centers = d["bin_centers"]
    bin_widths = d["bin_widths"]
    N = data["N"]

    logfuncdensity(function (params) 
        function bin_log_likelihood(i)
            xg = f(
                bin_centers[i], params.λ_g1, params.λ_g2, params.K_g, params.θ
            )
            expected_counts = bin_widths[i] * xg * N
            logpdf(Poisson(expected_counts), observed_counts[i])
        end

        idxs = eachindex(observed_counts)
        ll_value = bin_log_likelihood(idxs[1])
        for i in idxs[2:end]
            ll_value += bin_log_likelihood(i)
        end

        return ll_value
    end)
end

# The prior and likelihood can be passed to `BAT.jl` via the
# `PosteriorDensity`. We can then sample this density using `bat_sample()`.

posterior = PosteriorDensity(likelihood, prior);
samples = bat_sample(
    posterior, 
    MCMCSampling(mcalg=MetropolisHastings(), nsteps=10^4, nchains=4)
).result;

# The `SampledDensity` gives a quick overview of the results.

SampledDensity(posterior, samples)

#  ### Visualise results
# We can (roughly) check how well the fit reconstructs the truth with a simple comparison.

x_grid = range(0, stop=1, length=50)
sub_samples = bat_sample(samples, OrderedResampling(nsamples=200)).result

plot()
for i in eachindex(sub_samples)
    s = sub_samples[i].v
    xg = [xgx(x, s.λ_g1, s.λ_g2, s.K_g, s.θ) for x in bin_centers]
    plot!(bin_centers, xg .* bin_widths * N, alpha=0.1, lw=3, 
        color="darkorange", label="",)
end

xg = [xgx(x, λ_g1, λ_g2, K_g, θ) for x in bin_centers]
plot!(bin_centers, xg .* bin_widths * N, alpha=0.7, 
    label="Expected", lw=3, color="red")

scatter!(bin_centers, observed_counts, lw=3, label="Observed", color="black")

plot!(xlabel="x")

# In the above plot, the red line represents the truth, and the set
# of fainter lines represent samples from the posterior. 
#
# We can also look at marginal distributions for different parameters...
plot(
    samples, :(λ_g1),
    mean = true, std = true,
    nbins = 50, 
)

# ## Full model with all components
#
# We can now extend this approach to the full 9 components,
# there is not so much documentation here, as it follows the above case.
#
#
# ### Forward model
# Here the parametrisation of each component in decreasing order of importance.
# We also use the helper functions provided by `PartonDensity` to keep things tidy.

using PartonDensity

#

pdf_params = PDFParameters(λ_u=0.7, K_u=4.0, λ_d=0.5, K_d=6.0, 
                           λ_g1=0.7, λ_g2=-0.4, K_g=6.0, λ_q=-0.5,
                           weights=[1, 0.5, 0.3, 0.2, 0.1, 0.1, 0.1])

# Sanity check

int_xtotx(pdf_params)

# Plot true model

plot_input_pdfs(pdf_params)

# Generate example data

bins = 0.0:0.05:1.0
bin_widths = bins[2:end] - bins[1:end-1]
bin_centers = (bins[1:end-1] + bins[2:end]) / 2

N = 1000
nbins = size(bin_centers)[1]

expected_counts = zeros(nbins)
observed_counts = zeros(Integer, nbins)
for i in 1:nbins
    xt = xtotx(bin_centers[i], pdf_params) * N
    expected_counts[i] = bin_widths[i] * xt 
    observed_counts[i] = rand(Poisson(expected_counts[i]))
end

# Plot data and expectation

plot(bin_centers, [xtotx(x, pdf_params) for x in bin_centers] .* bin_widths * N,
     alpha=0.7, label="Expected", lw=3, color="red")
scatter!(bin_centers, observed_counts, lw=3, label="Observed", color="black")

# Store the data 

data = Dict()
data["N"] = N
data["bin_centers"] = bin_centers;
data["observed_counts"] = observed_counts;
data["bin_widths"] = bin_widths;

# ### Fit
#
# Prior

prior = NamedTupleDist(
    θ = Dirichlet(pdf_params.weights),
    λ_u = Truncated(Normal(pdf_params.λ_u, 0.5), 0, 1), #  Uniform(0, 1),
    K_u = Truncated(Normal(pdf_params.K_u, 1), 2, 10),
    λ_d = Truncated(Normal(pdf_params.λ_d, 0.5), 0, 1), # Uniform(0, 1),
    K_d = Truncated(Normal(pdf_params.K_d, 1), 2, 10),
    λ_g1 = Truncated(Normal(pdf_params.λ_g1, 1), 0, 1), 
    λ_g2 = Truncated(Normal(pdf_params.λ_g2, 1), -1, 0), 
    K_g = Truncated(Normal(pdf_params.K_g, 1), 2, 10),
    λ_q = Truncated(Normal(pdf_params.λ_q, 0.1), -1, 0),
);

# Likelihood

likelihood = let d = data, f = xtotx

    observed_counts = d["observed_counts"]
    bin_centers = d["bin_centers"]
    bin_widths = d["bin_widths"]
    N = data["N"]

    logfuncdensity(function (params) 
        function bin_log_likelihood(i)
            xt = f(bin_centers[i], params.λ_u, params.K_u, params.λ_d, params.K_d, 
                    params.λ_g1, params.λ_g2, params.K_g, params.λ_q, Vector(params.θ))
            expected_counts = bin_widths[i] * xt * N
            if expected_counts < 0
                expected_counts = 1e-3
            end
            logpdf(Poisson(expected_counts), observed_counts[i])
        end

        idxs = eachindex(observed_counts)
        ll_value = bin_log_likelihood(idxs[1])
        for i in idxs[2:end]
            ll_value += bin_log_likelihood(i)
        end

        return ll_value
    end)
end

# Run fit
#
# The next steps are commented for now as this hangs for some
# reason in the doc builds...

#posterior = PosteriorDensity(likelihood, prior);
#samples = bat_sample(posterior, MCMCSampling(mcalg=MetropolisHastings(), nsteps=10^4, nchains=2)).result;

# ### Visualise results

#x_grid = range(0, stop=1, length=50)
#sub_samples = bat_sample(samples, OrderedResampling(nsamples=200)).result

#plot()
"""
for i in eachindex(sub_samples)
    s = sub_samples[i].v
    xt = [xtotx(x, s.λ_u, s.K_u, s.λ_d, s.K_d, 
            s.λ_g1, s.λ_g2, s.K_g, s.λ_q, Vector(s.θ)) for x in bin_centers]
    plot!(bin_centers, xt .* bin_widths * N, alpha=0.1, lw=3, 
        color="darkorange", label="")
end
"""
#xt = [xtotx(x, pdf_params) for x in bin_centers]
#plot!(bin_centers, xt .* bin_widths * N, alpha=0.7, label="Expected", lw=3, color="red")

#scatter!(bin_centers, observed_counts, lw=3, label="Observed", color="black")
#plot!(xlabel="x")
nothing

# These first results are promising. We can also try changing the input
# parameters and priors to explore the performance of the fit.
