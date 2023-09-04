#module bla
using BAT, DensityInterface
using PartonDensity
using QCDNUM
using Plots, Colors , Random, Distributions, ValueShapes, ParallelProcessingTools
using StatsBase, LinearAlgebra
using DelimitedFiles
using ArgParse
export get_priors

function get_priors( parsed_args::Dict{String,Any} )
if parsed_args["parametrisation"] == "Dirichlet"
if (parsed_args["priorshift"]==0)
 #   println("seting prior from Shifted Prior set ",seedtxt)

prior = NamedTupleDist(
    θ = Dirichlet([20, 10, 20, 20, 5, 2.5, 1.5, 1.5, 0.5]),
    K_u = Truncated(Normal(3.5, 0.5), 1.0, 6.5),
    K_d = Truncated(Normal(3.5, 0.5), 1.0, 6.5),
    λ_g1 = Uniform(0., 1.),
    λ_g2 = Uniform(-1.0, -0.1),
    K_g =  Truncated(Normal(4., 1.5), 1.0, 8.5),
    λ_q = Uniform(-1.0, -0.1),
    K_q = Truncated(Normal(4., 1.5), 1.0, 9.5),
    Beta1 =  Truncated(Normal(0, 1), -5, 5),
    Beta2 =  Truncated(Normal(0, 1), -5, 5),
    beta0_1=  Truncated(Normal(0, 1), -5, 5), 
    beta0_2=   Truncated(Normal(0, 1), -5, 5),    
    beta0_3= Truncated(Normal(0, 1), -5, 5), 
    beta0_4=  Truncated(Normal(0, 1), -5, 5), 
    beta0_5=  Truncated(Normal(0, 1), -5, 5), 
    beta0_6=  Truncated(Normal(0, 1), -5, 5), 
    beta0_7=  Truncated(Normal(0, 1), -5, 5), 
    beta0_8=   Truncated(Normal(0, 1), -5, 5)
);
end

if (parsed_args["priorshift"]==1)
prior = NamedTupleDist(
    θ = Dirichlet([30, 10, 10, 10, 5, 2.5, 1.5, 1.5, 0.5]),
    K_u = Truncated(Normal(4.5, 0.5), 1.0, 6.5),
    K_d = Truncated(Normal(3.5, 0.5), 1.0, 6.5),
    λ_g1 = Uniform(0., 1.),
    λ_g2 = Uniform(-1.0, -0.1),
    K_g =  Truncated(Normal(4., 1.5), 1.0, 8.5),
    λ_q = Uniform(-1.0, -0.1),
    K_q = Truncated(Normal(4., 1.5), 1.0, 9.5),
    Beta1 =  Truncated(Normal(0, 1), -5, 5),
    Beta2 =  Truncated(Normal(0, 1), -5, 5),
    beta0_1=  Truncated(Normal(0, 1), -5, 5), 
    beta0_2=   Truncated(Normal(0, 1), -5, 5),    
    beta0_3= Truncated(Normal(0, 1), -5, 5), 
    beta0_4=  Truncated(Normal(0, 1), -5, 5), 
    beta0_5=  Truncated(Normal(0, 1), -5, 5), 
    beta0_6=  Truncated(Normal(0, 1), -5, 5), 
    beta0_7=  Truncated(Normal(0, 1), -5, 5), 
    beta0_8=   Truncated(Normal(0, 1), -5, 5)
);
elseif (parsed_args["priorshift"]==2)
  #  println("seting prior from Shifted Prior set ",seedtxt)
prior = NamedTupleDist(
    θ = Dirichlet([20, 10, 30, 30, 5, 2.5, 1.5, 1.5, 0.5]),
    K_u = Truncated(Normal(2.5, 0.5), 1.0, 6.5),
    K_d = Truncated(Normal(3.5, 0.5), 1.0, 6.5),
    λ_g1 = Uniform(0., 1.),
    λ_g2 = Uniform(-1.0, -0.1),
    K_g =  Truncated(Normal(4., 1.5),  1.0, 8.5),
    λ_q = Uniform(-1.0, -0.1),
    K_q = Truncated(Normal(4., 1.5), 1.0, 9.5),
    Beta1 =  Truncated(Normal(0, 1), -5, 5),
    Beta2 =  Truncated(Normal(0, 1), -5, 5),
    beta0_1=  Truncated(Normal(0, 1), -5, 5), 
    beta0_2=   Truncated(Normal(0, 1), -5, 5),    
    beta0_3= Truncated(Normal(0, 1), -5, 5), 
    beta0_4=  Truncated(Normal(0, 1), -5, 5), 
    beta0_5=  Truncated(Normal(0, 1), -5, 5), 
    beta0_6=  Truncated(Normal(0, 1), -5, 5), 
    beta0_7=  Truncated(Normal(0, 1), -5, 5), 
    beta0_8=   Truncated(Normal(0, 1), -5, 5)
);
end

end


if parsed_args["parametrisation"] == "Valence"
##FIXME!!!
weights = [5.0, 5.0, 1.0, 1.0, 1.0, 0.5, 0.5]
λ_u = 0.64;
K_u = 3.38;
λ_d = 0.67;
K_d = 4.73;
θ = get_θ_val( λ_u, K_u, λ_d, K_d, weights)
pdf_params = ValencePDFParams(λ_u=λ_u, K_u=K_u, λ_d=λ_d, K_d=K_d, λ_g1=0.50, λ_g2=-0.63, K_g=4.23, λ_q=-0.23, K_q=5.0, θ=θ);
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
    K_q=Truncated(Normal(pdf_params.K_q, 0.5), 3, 7),
    Beta1 =  Truncated(Normal(0, 1), -5, 5),
    Beta2 =  Truncated(Normal(0, 1), -5, 5),
    beta0_1=  Truncated(Normal(0, 1), -5, 5), 
    beta0_2=   Truncated(Normal(0, 1), -5, 5),    
    beta0_3= Truncated(Normal(0, 1), -5, 5), 
    beta0_4=  Truncated(Normal(0, 1), -5, 5), 
    beta0_5=  Truncated(Normal(0, 1), -5, 5), 
    beta0_6=  Truncated(Normal(0, 1), -5, 5), 
    beta0_7=  Truncated(Normal(0, 1), -5, 5), 
    beta0_8=   Truncated(Normal(0, 1), -5, 5)    
);
end

if parsed_args["parametrisation"] == "Bernstein"

prior = NamedTupleDist(
    θ = Dirichlet([28.0, 12.5, 20.0, 20.0, 10.0, 1.4, 0.2, 10.e-5, 0.3]),
#    initial_U = [Uniform(-10., 1.)],
    initial_U = [Truncated(Normal(-5, 3), -12, 5)],
    initial_D = [Uniform(10., 30.)],
    λ_g1 = Uniform(1., 2.0),
    λ_g2 = Uniform(-0.5, -0.3),
    K_g =  Uniform(5.,9.),
    λ_q = Uniform(-0.5, -0.),
    K_q = Uniform(3., 7.),
    bspoly_params = [0,5,1,5,0,6],
    Beta1 =  Truncated(Normal(0, 1), -5, 5),
    Beta2 =  Truncated(Normal(0, 1), -5, 5),
    beta0_1=  Truncated(Normal(0, 1), -5, 5), 
    beta0_2=   Truncated(Normal(0, 1), -5, 5),    
    beta0_3= Truncated(Normal(0, 1), -5, 5), 
    beta0_4=  Truncated(Normal(0, 1), -5, 5), 
    beta0_5=  Truncated(Normal(0, 1), -5, 5), 
    beta0_6=  Truncated(Normal(0, 1), -5, 5), 
    beta0_7=  Truncated(Normal(0, 1), -5, 5), 
    beta0_8=   Truncated(Normal(0, 1), -5, 5)
    )
end

return prior

end

#end
