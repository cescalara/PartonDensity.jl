using Test

@testset "PartonDensity" begin

    include("test_parametrisations.jl")
    include("test_forward_model.jl")
    
end
