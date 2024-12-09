import PartonDensity

@testset "genpseudodata and examples" begin
    function run_with_args(jl_filename, args)
                pkgdir = dirname(dirname(pathof(PartonDensity)))
                old_args = collect(ARGS)
                empty!(ARGS)
                append!(ARGS, args)
                include(joinpath(pkgdir, "utils", jl_filename))
                empty!(ARGS)
                append!(ARGS, old_args)
    end

    mktempdir(prefix = "test_QCDNUM_jl") do dir
        cd(dir) do
            mkdir("pseudodata")
            run_with_args("generatepseudodata.jl", ["-s", "42", "-p", "Dirichlet", "-f", "0.5", "-c", "true"])
            run_with_args("generatepseudodata.jl", ["-s", "43", "-p", "Dirichlet", "-f", "1.5", "-c", "true"])

            mkdir("fitresults")
            run_with_args("PDFfit.jl", ["-s", "45", "-p", "Dirichlet", "-d", "simulation-Dirichlet-42", "simulation-Dirichlet-43", "-n", "250", "--max_ncycles=0", "--nsteps_per_cycle=10", "--nsteps_final=10", "--strict=false", "--dummylikelihood=true"])

            mkdir("CABCHSV2023nnlo")
            run_with_args("generatelhapdf.jl", ["-s", "43", "-p", "Dirichlet", "-f", "fit-Dirichlet-0-45-simulation-Dirichlet-42simulation-Dirichlet-43"])
        end
    end
end
