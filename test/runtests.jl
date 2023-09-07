using Test
using ArgCheck
include("../data/ZEUS_I1787035/ZEUS_I1787035.jl")
@testset "PartonDensity" begin

    include("test_parametrisations.jl")
    include("test_forward_model.jl")
    include("test_posterior.jl")
end
