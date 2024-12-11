# This file is a part of PartonDensity.jl, licensed under the MIT License (MIT).

module PartonDensityBATExt

using PartonDensity

using BAT

using DensityInterface
using Distributions
using ValueShapes
using QCDNUM

include("BAT/fit.jl")
include("BAT/bernstein_fit.jl")

end # module PartonDensityBATExt
