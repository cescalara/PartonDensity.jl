using Test
using ArgCheck
const MD_ZEUS_I1787035=include("../data/ZEUS_I1787035/ZEUS_I1787035.jl")
@testset "PartonDensity" begin

    include("test_parametrisations.jl")
    include("test_forward_model.jl")
    include("test_posterior.jl")
end
