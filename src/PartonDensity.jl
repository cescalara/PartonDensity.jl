module PartonDensity

include("parametrisation.jl")
include("qcdnum_interface.jl")
include("cross_section.jl")
include("forward_model.jl")
include("forward_model_sysErr.jl")
include("zeus.jl")
include("Bernsteinparametrisation.jl")

using Requires
function __init__()
    @require BAT = "c0cd4b16-88b7-57fa-983b-ab80aecada7e" include("fit.jl")
end

end
