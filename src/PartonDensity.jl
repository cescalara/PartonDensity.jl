module PartonDensity

using ArgCheck
using Distributions
using DocStringExtensions
using Random
using Parameters
using Plots

include("parton.jl")
include("parametrisations/parametrisations.jl")
include("qcdnum_interface.jl")
include("cross_section.jl")
include("forward_model.jl")
include("bernstein_forward_model.jl")

using Requires
function __init__()
    @require BAT = "c0cd4b16-88b7-57fa-983b-ab80aecada7e" include("fit.jl")
    @require BAT = "c0cd4b16-88b7-57fa-983b-ab80aecada7e" include("bernstein_fit.jl")
end

end