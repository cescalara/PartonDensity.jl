#!/usr/bin/julia
using BAT, DensityInterface
using PartonDensity
using QCDNUM
using Random, Distributions, ValueShapes, ParallelProcessingTools
using StatsBase, LinearAlgebra
using DelimitedFiles
using ArgParse
import HDF5
 
include("priors.jl")
include(string(dirname(pathof(PartonDensity)),"/../data/ZEUS_I1787035/ZEUS_I1787035.jl"))

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
            
        "--max_ncycles"
            help = "max_ncycles"
            arg_type = Int
            default = 50
            
         "--nsteps_per_cycle"
            help = "nsteps_per_cycle"
            arg_type = Int
            default = 10000
         "--nsteps_final"
            help = "nsteps_final"
            arg_type = Int
            default = 1000

         "--strict"
            help = "nsteps_final"
            arg_type = Bool
            default = true

         "--dummylikelihood"
            help = "dummylikelihood"
            arg_type = Bool
            default = false

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

prior=get_priors(parsed_args)


sim_data = Dict{String,Any}()

sim_data["counts_obs_ep"]=MD_G.m_Data_Events_ePp
sim_data["counts_obs_em"]=MD_G.m_Data_Events_eMp
MD_TEMP::MetaData = MD_G

if parsed_args["pseudodata"] != "data"
  somepdf_params, sim_data, MD_TEMP = pd_read_sim(string("pseudodata/",parsed_args["pseudodata"],".h5"),MD_G);
end
MD_LOCAL::MetaData = MD_TEMP



likelihood = let d = sim_data
counts_obs_ep = d["counts_obs_ep"]
counts_obs_em = d["counts_obs_em"]
nbins = size(d["counts_obs_ep"])[1]
logfuncdensity(function (params)
       if parsed_args["dummylikelihood"]
         return -100.0;
       end
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
         counts_pred_ep, counts_pred_em = @critical  forward_model(pdf_params, qcdnum_params, splint_params, quark_coeffs,MD_LOCAL,ParErrs );
       end

       if parsed_args["parametrisation"] == "Dirichlet"
         pdf_params = DirichletPDFParams(K_u=params.K_u, K_d=params.K_d, λ_g1=params.λ_g1, λ_g2=params.λ_g2, K_g=params.K_g, λ_q=params.λ_q, K_q=params.K_q, θ=params.θ)
         ParErrs = [params.beta0_1,params.beta0_2,params.beta0_3,params.beta0_4, params.beta0_5,params.beta0_6,params.beta0_7,params.beta0_8]
         counts_pred_ep, counts_pred_em = @critical  forward_model(pdf_params, qcdnum_params, splint_params, quark_coeffs,MD_LOCAL,ParErrs );
       end

       if parsed_args["parametrisation"] == "Valence"
         θ = get_scaled_θ(params.λ_u, params.K_u, params.λ_d,params.K_d, Vector(params.θ_tmp))
         pdf_params = ValencePDFParams(λ_u=params.λ_u, K_u=params.K_u, λ_d=params.λ_d, K_d=params.K_d, λ_g1=params.λ_g1, λ_g2=params.λ_g2, K_g=params.K_g, λ_q=params.λ_q, K_q=params.K_q, θ=θ)
         ParErrs = [params.beta0_1,params.beta0_2,params.beta0_3,params.beta0_4, params.beta0_5,params.beta0_6,params.beta0_7,params.beta0_8]
         counts_pred_ep, counts_pred_em = @critical  forward_model(pdf_params, qcdnum_params, splint_params, quark_coeffs,MD_LOCAL,ParErrs );
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
                counts_pred_ep[i] =counts_pred_ep[i]*(1+MD_LOCAL.Ld_ePp_uncertainty*params.Beta1)
                counts_pred_em[i] =counts_pred_em[i]*(1+MD_LOCAL.Ld_eMp_uncertainty*params.Beta2)                
                ll_value += logpdf(Poisson(counts_pred_ep[i]), counts_obs_ep[i])
                ll_value += logpdf(Poisson(counts_pred_em[i]), counts_obs_em[i])
            end
            return ll_value
    end)
end


posterior = PosteriorDensity(likelihood, prior);
mcalg = MetropolisHastings(proposal=BAT.MvTDistProposal(10.0))
convergence = BrooksGelmanConvergence(threshold=1.3);
burnin = MCMCMultiCycleBurnin(max_ncycles=parsed_args["max_ncycles"],nsteps_per_cycle=parsed_args["nsteps_per_cycle"],nsteps_final=parsed_args["nsteps_final"]);
samples = bat_sample(posterior, MCMCSampling(mcalg=mcalg, nsteps=parsed_args["nsteps"], nchains=parsed_args["nchains"],strict=parsed_args["strict"])).result;

fname=string("fitresults/fit-",parsed_args["parametrisation"],"-",parsed_args["priorshift"],"-",seedtxt,"-",parsed_args["pseudodata"])
bat_write(string(fname,".h5"), samples)
QCDNUM.save_params(string(fname,"_qcdnum.h5"), qcdnum_params)
end 

main()



