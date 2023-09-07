using Test
using ArgCheck

@testset "PartonDensity" begin

    include("test_parametrisations.jl")
    include("test_forward_model.jl")
    include("test_posterior.jl")
end
