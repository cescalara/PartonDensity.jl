#!/usr/bin/julia
using BAT, DensityInterface
using PartonDensity
using QCDNUM
using Plots, Colors , Random, Distributions, ValueShapes, ParallelProcessingTools
using StatsBase, LinearAlgebra

using ArgParse

function parse_commandline()
    s = ArgParseSettings()

    @add_arg_table s begin
        "--seed", "-s"
            help = "Seed"
            arg_type = Int
            default = 42
        "--parametrisation", "-p"
            help = "Parametrisation -- Dirichlet or Valence"
            arg_type = String
            default = "Dirichlet"
        "--dataset", "-d"
            help = "Dataset ID to generate the data"
            arg_type = String
            default = "../data/ZEUS_I1787035/ZEUS_I1787035.jl"
        "--lumifactor", "-f"
            help = "Lumi factor for the Dataset"
            arg_type = Float64
            default = 1.0

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

if parsed_args["parametrisation"] == "Dirichlet"
  θ = [ 0.228, 0.104, 0.249, 0.249, 0.104, 0.052, 0.010, 0.005, 0.0005]
  θ_sum=sum(θ[1:9])
  θ=θ/θ_sum  
  pdf_params = DirichletPDFParams(K_u=3.7, K_d=3.7, λ_g1=0.5, λ_g2=-0.5, K_g=5.0,λ_q=-0.5, K_q=6.0, θ=θ);
end
if parsed_args["parametrisation"] == "Valence"
#  weights = [5.0, 5.0, 1.0, 1.0, 1.0, 0.5, 0.5]
#  θ = get_θ_val(rng, λ_u, K_u, λ_d, K_d, weights)
  λ_u = 0.64;
  K_u = 3.38;
  λ_d = 0.67;
  K_d = 4.73;
  θ=[0.22, 0.10, 0.24, 0.24, 0.10,0.05, 0.01, 0.005, 0.0005]
  θ_sum=sum(θ[1:9])
  θ=θ/θ_sum
  pdf_params = ValencePDFParams(λ_u=λ_u, K_u=K_u, λ_d=λ_d, K_d=K_d, λ_g1=0.50, λ_g2=-0.63, K_g=4.23, λ_q=-0.23, K_q=5.0, θ=θ);
end


if parsed_args["parametrisation"] == "Bernstein"
  θ=[.33, .13, .27, .17, .073, 0.014, 0.002, 0.000001, .003]
  θ_sum=sum(θ[1:9])
  θ=θ/θ_sum  
  bspoly_params = [1,4,0,4,0,5];
  λ_g1=1.5;
  λ_g2=-0.4;
  K_g=6.0;
  λ_q=-0.25;
  K_q=5.0;
  initial_U = [-8.];
  initial_D = [15.0];
  vec_bspp = Vector(bspoly_params)
  bspoly_params = [[vec_bspp[Int(2 * i - 1)], vec_bspp[Int(2 * i)]] for i in 1:length(vec_bspp)/2]
  bspoly_params_d = 0
  try
    vec_bsppd = Vector(bspoly_params_d)
    bspoly_params_d = [[vec_bsppd[Int(2 * i - 1)], vec_bsppd[Int(2 * i)]] for i in 1:length(vec_bsppd)/2]
  catch err
    bspoly_params_d = bspoly_params
  end 
  initU = Vector(initial_U)
  initD = Vector(initial_D)       
  pdf_params = BernsteinDirichletPDFParams(initial_U=initU,
                    initial_D=initD,
                    λ_g1=λ_g1, λ_g2=λ_g2,
                    K_g=K_g, λ_q=λ_q, K_q=K_q,
                    bspoly_params=bspoly_params,
                    bspoly_params_d=bspoly_params_d,
                    θ=Vector(θ)
  )
end


qcdnum_grid = QCDNUM.GridParams(x_min=[1.0e-3, 1.0e-1, 5.0e-1], x_weights=[1, 2, 2], nx=100, qq_bounds=[1.0e2, 3.0e4], qq_weights=[1.0, 1.0], nq=50, spline_interp=3)
qcdnum_params = QCDNUM.EvolutionParams(order=2, α_S=0.118, q0=100.0, grid_params=qcdnum_grid,n_fixed_flav=5, iqc=1, iqb=1, iqt=1, weight_type=1);
    
splint_params = QCDNUM.SPLINTParams();
quark_coeffs = QuarkCoefficients();


include(parsed_args["dataset"])
MD_LOCAL::MetaData =  MetaData(MD_G.name, 
                               MD_G.Ld_ePp*parsed_args["lumifactor"] , 
                               MD_G.Ld_eMp*parsed_args["lumifactor"], 
                               MD_G.Ld_ePp_uncertainty,
                               MD_G.Ld_eMp_uncertainty, 
                               MD_G.sqrtS,
                               MD_G.nsyst,
                               PD_DUMMY)


forward_model_init(qcdnum_params, splint_params)
counts_pred_ep, counts_pred_em = forward_model(pdf_params, qcdnum_params,splint_params, quark_coeffs,MD_LOCAL);
    
nbins = size(counts_pred_ep)[1]
counts_obs_ep = zeros(UInt64, nbins)
counts_obs_em = zeros(UInt64, nbins)

for i in 1:nbins
    counts_obs_ep[i] = rand(rng,Poisson(counts_pred_ep[i]))
    counts_obs_em[i] = rand(rng,Poisson(counts_pred_em[i]))
end

sim_data = Dict{String,Any}()
sim_data["nbins"] = nbins;
sim_data["counts_obs_ep"] = counts_obs_ep;
sim_data["counts_obs_em"] = counts_obs_em;

pd_write_sim(string("pseudodata/simulation-",parsed_args["parametrisation"],"-",seedtxt,".h5"), pdf_params, sim_data, MD_G)
end

main()
