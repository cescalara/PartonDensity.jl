push!(LOAD_PATH,"../src/")

using Documenter, PartonDensity

makedocs(modules=[PartonDensity], sitename="PartonDensity.jl")

deploydocs(
    devbranch = "main",
    repo = "github.com/cescalara/PartonDensity.jl.git",
)
