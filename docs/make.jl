push!(LOAD_PATH,"../src/")

using Documenter, PartonDensity
using Literate

# Generate examples notebooks

gen_content_dir = joinpath(@__DIR__, "src")

pdf_parametrisation_src = joinpath(@__DIR__, "..", "examples", "pdf_parametrisation.jl")
forward_model_src = joinpath(@__DIR__, "..", "examples", "forward_model.jl")
transfer_matrix_src = joinpath(@__DIR__, "..", "examples", "transfer_matrix.jl")

Literate.markdown(pdf_parametrisation_src, gen_content_dir, name="pdf_parametrisation")
Literate.markdown(forward_model_src, gen_content_dir, name="forward_model")
Literate.markdown(transfer_matrix_src, gen_content_dir, name="transfer_matrix")

Introduction = "Introduction" => "index.md"

Examples = "Examples" => ["pdf_parametrisation.md",
                          "forward_model.md", "transfer_matrix.md"]

API = "API" => "api.md"

PAGES = [Introduction, Examples, API]

makedocs(modules=[PartonDensity], sitename="PartonDensity.jl", pages=PAGES)

deploydocs(
    devbranch = "main",
    repo = "github.com/cescalara/PartonDensity.jl.git",
)
