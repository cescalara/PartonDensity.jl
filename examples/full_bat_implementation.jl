# # A fit with BAT.jl

using Pkg
Pkg.activate("../docs")
#Pkg.add(url="https://github.com/bat/BAT.jl.git")
#Pkg.add(url="https://github.com/cescalara/QCDNUM.jl.git", rev="cfuncs")
#Pkg.develop(PackageSpec(path="../")) 
#Pkg.instantiate()

using BAT, DensityInterface
using PartonDensity
using QCDNUM
using Plots, Random, Distributions, ValueShapes, ParallelProcessingTools

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
# definition of the input PDFs.

pdf_params = PDFParameters(λ_u=0.7, K_u=4.0, λ_d=0.5, K_d=6.0,
    λ_g1=-0.4, λ_g2=0.7, K_g=6.0, λ_q=0.5, weights=[1, 0.5, 0.3, 0.2, 0.1, 0.1, 0.1]);

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
splint_params = SPLINTParameters(nuser=1000);
quark_coeffs = QuarkCoefficients();

# initialise QCDNUM
forward_model_init(qcdnum_grid, qcdnum_params)

# run forward model 
counts_pred_ep, counts_pred_em = forward_model(pdf_params, qcdnum_params, 
    splint_params, quark_coeffs);

# ## Debugging...
#
# The implementation below crashes. It seems even calling forward_model in a loop for >> 10 times has similar issues. 
#
# Ideas
# * Memory issue? (Test reinitialising QCDNUM on each iteration...)
# * Rare problem - need to debug with right tools

for i in 1:100
    #forward_model_init(qcdnum_grid, qcdnum_params)
    counts_pred_ep, counts_pred_em = forward_model(pdf_params, qcdnum_params, 
    splint_params, quark_coeffs);
end

# ## Continue simulating data...
#
# take a poisson sample
# nbins = size(counts_pred_ep)[1]
# counts_obs_ep = zeros(Integer, nbins)
# counts_obs_em = zeros(Integer, nbins)

# for i in 1:nbins
#     counts_obs_ep[i] = rand(Poisson(counts_pred_ep[i]))
#     counts_obs_em[i] = rand(Poisson(counts_pred_em[i]))
# end

# #

# plot(1:nbins, counts_pred_ep, label="Expected counts (eP)", color="blue")
# plot!(1:nbins, counts_pred_em, label="Expected counts (eM)", color="red")
# scatter!(1:nbins, counts_obs_ep, label="Detected counts (eP)", color="blue")
# scatter!(1:nbins, counts_obs_em, label="Detected counts (eM)", color="red")
# plot!(xlabel="Bin number")

# # store
# sim_data = Dict()
# sim_data["nbins"] = nbins;
# sim_data["counts_obs_ep"] = counts_obs_ep;
# sim_data["counts_obs_em"] = counts_obs_em;

# ## Fit the simulated data
#
# Now we can try to fit this simulated data using `Bat.jl`.
# The first step is to define the prior and likelihood.

# prior = NamedTupleDist(
#     θ = Dirichlet([1, 0.5, 0.3, 0.2, 0.1, 0.1, 0.1]),
#     λ_u = Truncated(Normal(pdf_params.λ_u, 0.05), 0, 1), 
#     K_u = Uniform(3.5, 4.5),
#     λ_d = Truncated(Normal(pdf_params.λ_d, 0.05), 0, 1), 
#     K_d = Uniform(5.5, 6.5),
#     λ_g1 = Uniform(-0.45, -0.35),
#     λ_g2 = Uniform(0.6, 0.8),
#     K_g =  Uniform(5.5, 6.5),
#     λ_q = Truncated(Normal(pdf_params.λ_q, 0.05), -1, 0),
# );

# #

# likelihood = let d = sim_data

#     counts_obs_ep = d["counts_obs_ep"]
#     counts_obs_em = d["counts_obs_em"]
#     nbins = d["nbins"]

#     logfuncdensity(function (params)
            
#             pdf_params = PDFParameters(λ_u=params.λ_u, K_u=params.K_u, λ_d=params.λ_d,
#                 K_d=params.K_d, λ_g1=params.λ_g1, λ_g2=params.λ_g2, K_g=params.K_g,
#                 λ_q=params.λ_q, θ=Vector(params.θ))
#             @critical counts_pred_ep, counts_pred_em = forward_model(pdf_params, 
#                 qcdnum_params, splint_params, quark_coeffs);

#             ll_value = 0.0
#             for i in 1:nbins
#                 ll_value += logpdf(Poisson(counts_pred_ep[i]), counts_obs_ep[i])
#                 ll_value += logpdf(Poisson(counts_pred_em[i]), counts_obs_em[i])
#             end

#             return ll_value
#     end)
# end

# #

# posterior = PosteriorDensity(likelihood, prior);
# samples = bat_sample(
#     posterior,
#     MCMCSampling(mcalg=MetropolisHastings(), nsteps=100, nchains=1)
# ).result;
