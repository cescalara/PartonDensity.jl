push!(LOAD_PATH,"../src/")

using Documenter, PartonDensity

Introduction = "Introduction" => "index.md"

Examples = "Examples" => "examples.md"

API = "API" => "api.md"

PAGES = [Introduction, Examples, API]

makedocs(modules=[PartonDensity], sitename="PartonDensity.jl", pages=PAGES)

deploydocs(
    devbranch = "main",
    repo = "github.com/cescalara/PartonDensity.jl.git",
)
