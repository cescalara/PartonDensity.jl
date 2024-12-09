push!(LOAD_PATH, "../src/")

using Documenter, PartonDensity
using Literate
import BAT

# Doctest setup
DocMeta.setdocmeta!(
    PartonDensity,
    :DocTestSetup,
    :(using PartonDensity);
    recursive=true,
)

# Generate examples notebooks

gen_content_dir = joinpath(@__DIR__, "src")

pdf_parametrisation_src = joinpath(@__DIR__, "..", "examples", "pdf_parametrisation.jl")
zeus_interface_src = joinpath(@__DIR__, "..", "examples", "zeus_interface.jl")
forward_model_src = joinpath(@__DIR__, "..", "examples", "forward_model.jl")
#bat_valence_src = joinpath(@__DIR__, "..", "examples", "bat_fit_valence.jl")
#bat_dirichlet_src = joinpath(@__DIR__, "..", "examples", "bat_fit_dirichlet.jl")
#bat_dirichlet_syserr_src = joinpath(@__DIR__, "..", "examples", "bat_fit_dirichlet_syserr.jl")

Literate.markdown(pdf_parametrisation_src, gen_content_dir, name="pdf_parametrisation")
Literate.markdown(zeus_interface_src, gen_content_dir, name="zeus_interface")
Literate.markdown(forward_model_src, gen_content_dir, name="forward_model")
#Literate.markdown(bat_valence_src, gen_content_dir, name="bat_fit_valence")
#Literate.markdown(bat_dirichlet_src, gen_content_dir, name="bat_fit_dirichlet")
#Literate.markdown(bat_dirichlet_syserr_src, gen_content_dir, name="bat_fit_dirichlet_syserr")

Examples = "Examples" => [
    "pdf_parametrisation.md",
    "zeus_interface.md", "forward_model.md",
    # "bat_fit_valence.md", "bat_fit_dirichlet.md", "bat_fit_dirichlet_syserr.md"
]

makedocs(
    sitename="PartonDensity.jl",
    modules=[PartonDensity],
    pages= [
        "Introduction" => "index.md",
        Examples,
        "API" => "api.md",
    ],
    format = Documenter.HTML(
        prettyurls = !("local" in ARGS),
        canonical = "https://cescalara.github.io/PartonDensity.jl/stable/"
    ),
    doctest = ("fixdoctests" in ARGS) ? :fix : true,
    linkcheck = !("nonstrict" in ARGS),
    warnonly = ("nonstrict" in ARGS),
)

deploydocs(
    repo="github.com/cescalara/PartonDensity.jl.git",
    forcepush = true,
    push_preview = true,
)
