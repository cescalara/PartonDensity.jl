# This file is a part of PartonDensity.jl, licensed under the MIT License (MIT).

import Test
import Aqua
import PartonDensity

Test.@testset "Package ambiguities" begin
    Test.@test isempty(Test.detect_ambiguities(PartonDensity))
end # testset

Test.@testset "Aqua tests" begin
    Aqua.test_all(
        PartonDensity,
        ambiguities = true
    )
end # testset
