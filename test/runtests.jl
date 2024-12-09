using Test

const MD_ZEUS_I1787035=include("../data/ZEUS_I1787035/ZEUS_I1787035.jl")

@testset "PartonDensity" begin
    include("test_aqua.jl")
    include("test_parametrisations.jl")
    include("test_forward_model.jl")
    include("test_posterior.jl")
    include("test_docs.jl")

    include("test_generatepseudodata.jl")
end
