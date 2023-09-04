#!/usr/bin/julia
using BAT, DensityInterface
using PartonDensity
using QCDNUM
using Plots, Colors , Random, Distributions, ValueShapes, ParallelProcessingTools
using StatsBase, LinearAlgebra
using DelimitedFiles
using ArgParse
import HDF5
 
include("priors.jl")


function parse_commandline()
    s = ArgParseSettings()

    @add_arg_table s begin
        "--priorshift"
            help = "Shift (variation) of priors. Only for Dirichlet"
            arg_type = Int
            default = 0
        "--seed", "-s"
            help = "Seed"
            arg_type = Int
            default = 42
        "--nchains", "-c"
            help = "Chains"
            arg_type = Int
            default = 2

        "--nsteps", "-n"
            help = "Number of steps"
            arg_type = Int
            default = 10^3
        "--parametrisation", "-p"
            help = "Parametrisation -- Dirichlet or Valence"
            arg_type = String
            default = "Dirichlet"
        "--pseudodata", "-d"
            help = "Input pseudodata -- file in the pseudodata directory w/o the extension"
            arg_type = String
            default = "data"
    end

    return parse_args(s)
end

function main()
    parsed_args = parse_commandline()
    println("Parsed args:")
    for (arg,val) in parsed_args
        println("  $arg  =>  $val")
    end

gr(fmt=:png);


counts_obs_ep = get_data_events(0)
counts_obs_em = get_data_events(1)

nbins = size(counts_obs_ep)[1]

seed=parsed_args["seed"]
println(seed)
seedtxt=string(seed)
rng = MersenneTwister(seed)


# first specify QCDNUM inputs
qcdnum_grid = QCDNUM.GridParams(x_min=[1.0e-3, 1.0e-1, 5.0e-1], x_weights=[1, 2, 2], nx=100, qq_bounds=[1.0e2, 3.0e4], qq_weights=[1.0, 1.0], nq=50, spline_interp=3)
qcdnum_params = QCDNUM.EvolutionParams(order=2, α_S=0.118, q0=100.0, grid_params=qcdnum_grid, n_fixed_flav=5, iqc=1, iqb=1, iqt=1, weight_type=1);


# now SPLINT and quark coefficients
splint_params = QCDNUM.SPLINTParams();
quark_coeffs = QuarkCoefficients();

# initialise QCDNUM
forward_model_init(qcdnum_params, splint_params)


sim_data = Dict{String,Any}()
sim_data["nbins"] = nbins;
sim_data["counts_obs_ep"]=get_data_events(0)
sim_data["counts_obs_em"]=get_data_events(1)
θ = [ 0.228, 0.104, 0.249, 0.249, 0.104, 0.052, 0.010, 0.005, 0.0005]
θ_sum=sum(θ[1:9])
θ=θ/θ_sum  
somepdf_params = DirichletPDFParams(K_u=3.7, K_d=3.7, λ_g1=0.5, λ_g2=-0.5, K_g=5.0,λ_q=-0.5, K_q=6.0, θ=θ);

if parsed_args["pseudodata"] != "data"
  somepdf_params, sim_data = pd_read_sim(string("pseudodata/",parsed_args["pseudodata"],".h5"));
end

prior=get_priors(parsed_args)

# The `@critical` macro is used because `forward_model()` is currently not thread safe, so
# this protects it from being run in parallel.


likelihood = let d = sim_data
counts_obs_ep = d["counts_obs_ep"]
counts_obs_em = d["counts_obs_em"]
nbins = d["nbins"]
logfuncdensity(function (params)
       if parsed_args["parametrisation"] == "Bernstein"
         vec_bspp = Vector(params.bspoly_params)
         bspoly_params = [[vec_bspp[Int(2 * i - 1)], vec_bspp[Int(2 * i)]] for i in 1:length(vec_bspp)/2]
         bspoly_params_d = 0
         try
           vec_bsppd = Vector(params.bspoly_params_d)
           bspoly_params_d = [[vec_bsppd[Int(2 * i - 1)], vec_bsppd[Int(2 * i)]] for i in 1:length(vec_bsppd)/2]
         catch err
                bspoly_params_d = bspoly_params
         end
         initU = Vector(params.initial_U)
         initD = Vector(params.initial_D)
         pdf_params = BernsteinDirichletPDFParams(initial_U=initU,initial_D=initD, λ_g1=params.λ_g1, λ_g2=params.λ_g2,
                    K_g=params.K_g, λ_q=params.λ_q, K_q=params.K_q,bspoly_params=bspoly_params,bspoly_params_d=bspoly_params_d,θ=Vector(params.θ))
         ParErrs = [params.beta0_1,params.beta0_2,params.beta0_3,params.beta0_4, params.beta0_5,params.beta0_6,params.beta0_7,params.beta0_8]
         counts_pred_ep, counts_pred_em = @critical  forward_model(pdf_params, qcdnum_params, splint_params, quark_coeffs,ParErrs );
       end

       if parsed_args["parametrisation"] == "Dirichlet"
         pdf_params = DirichletPDFParams(K_u=params.K_u, K_d=params.K_d, λ_g1=params.λ_g1, λ_g2=params.λ_g2, K_g=params.K_g, λ_q=params.λ_q, K_q=params.K_q, θ=params.θ)
         ParErrs = [params.beta0_1,params.beta0_2,params.beta0_3,params.beta0_4, params.beta0_5,params.beta0_6,params.beta0_7,params.beta0_8]
         counts_pred_ep, counts_pred_em = @critical  forward_model(pdf_params, qcdnum_params, splint_params, quark_coeffs,ParErrs );
       end

       if parsed_args["parametrisation"] == "Valence"
         θ = get_scaled_θ(params.λ_u, params.K_u, params.λ_d,params.K_d, Vector(params.θ_tmp))
         pdf_params = ValencePDFParams(λ_u=params.λ_u, K_u=params.K_u, λ_d=params.λ_d, K_d=params.K_d, λ_g1=params.λ_g1, λ_g2=params.λ_g2, K_g=params.K_g, λ_q=params.λ_q, K_q=params.K_q, θ=θ)
         ParErrs = [params.beta0_1,params.beta0_2,params.beta0_3,params.beta0_4, params.beta0_5,params.beta0_6,params.beta0_7,params.beta0_8]
         counts_pred_ep, counts_pred_em = @critical  forward_model(pdf_params, qcdnum_params, splint_params, quark_coeffs,ParErrs );
       end
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
                counts_pred_ep[i] =counts_pred_ep[i]*(1+0.018*params.Beta1)
                counts_pred_em[i] =counts_pred_em[i]*(1+0.018*params.Beta2)                
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
mcalg = MetropolisHastings(proposal=BAT.MvTDistProposal(10.0))
convergence = BrooksGelmanConvergence(threshold=1.3);
burnin = MCMCMultiCycleBurnin(max_ncycles=50);
samples = bat_sample(posterior, MCMCSampling(mcalg=mcalg, nsteps=parsed_args["nsteps"], nchains=parsed_args["nchains"])).result;
# Let's save the result for further analysis
bat_write(string("fitresults/fit-",parsed_args["parametrisation"],"-",parsed_args["priorshift"],"-",seedtxt,"-",parsed_args["pseudodata"],".h5"), samples)
end 

main()



