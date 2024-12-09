# This file is a part of PartonDensity.jl, licensed under the MIT License (MIT).

using Test
using PartonDensity
import Documenter

Documenter.DocMeta.setdocmeta!(
    PartonDensity,
    :DocTestSetup,
    :(using PartonDensity);
    recursive=true,
)
Documenter.doctest(PartonDensity)
