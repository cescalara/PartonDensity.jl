module PartonDensity

using ArgCheck
using Distributions
using DocStringExtensions
using Random
using Parameters

include("parton.jl")
include("parametrisations/parametrisations.jl")
include("metadata.jl")
include("parameters_IO.jl")
include("cross_section.jl")
include("forward_model.jl")

include("extdefs_BAT.jl")
include("extdefs_Plots.jl")

include("register_extdeps.jl")

function __init__()
    _register_extension_deps(
        get_prior => :BAT,
        get_likelihood => :BAT,
        plot_input_pdfs => :Plots,
        plot_model_space => :Plots,
        plot_data_space => :Plots,
    )
end

end
