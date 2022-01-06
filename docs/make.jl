push!(LOAD_PATH,"../src/")

using Documenter, PartonDensity
using Literate

# Generate examples notebooks

gen_content_dir = joinpath(@__DIR__, "src")

pdf_parametrisation_src = joinpath(@__DIR__, "..", "examples", "pdf_parametrisation.jl")
transfer_matrix_src = joinpath(@__DIR__, "..", "examples", "transfer_matrix.jl")
forward_model_src = joinpath(@__DIR__, "..", "examples", "forward_model.jl")
prior_bat_src = joinpath(@__DIR__, "..", "examples", "prior_bat_implementation.jl")

Literate.markdown(pdf_parametrisation_src, gen_content_dir, name="pdf_parametrisation")
Literate.markdown(transfer_matrix_src, gen_content_dir, name="transfer_matrix")
Literate.markdown(forward_model_src, gen_content_dir, name="forward_model")
Literate.markdown(prior_bat_src, gen_content_dir, name="prior_bat_implementation")

Introduction = "Introduction" => "index.md"

Examples = "Examples" => ["pdf_parametrisation.md",
                          "transfer_matrix.md", "forward_model.md",
                          "prior_bat_implementation.md"]

API = "API" => "api.md"

PAGES = [Introduction, Examples, API]

makedocs(modules=[PartonDensity], sitename="PartonDensity.jl", pages=PAGES)

deploydocs(
    devbranch = "main",
    repo = "github.com/cescalara/PartonDensity.jl.git",
)
