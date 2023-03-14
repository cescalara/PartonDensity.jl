using PartonDensity
using Test
using Random, Distributions
using Plots, ValueShapes
import BAT

@testset "Valence PDF parametrisation" begin

    local val_pdf_params

    for i in 1:100

        λ_u = rand(Uniform(0, 1))
        K_u = rand(Uniform(2, 10))
        λ_d = rand(Uniform(0, 1))
        K_d = rand(Uniform(2, 10))
        λ_g1 = rand(Uniform(0, 1))
        λ_g2 = rand(Uniform(-1, 0))
        K_g = rand(Uniform(2, 10))
        λ_q = rand(Uniform(-1, 0))
        K_q = rand(Uniform(1, 5))
        weights = ones(7)
        θ = get_θ_val(λ_u, K_u, λ_d, K_d, weights)

        val_pdf_params = ValencePDFParams(λ_u=λ_u, K_u=K_u,
            λ_d=λ_d, K_d=K_d, λ_g1=λ_g2,
            λ_g2=λ_g2, K_g=K_g, λ_q=λ_q, K_q=K_q,
            θ=θ)

        @test int_xtotx(val_pdf_params) ≈ 1.0

    end

    @test typeof(val_pdf_params) == ValencePDFParams{Float64,Vector{Float64}}

    p = plot_input_pdfs(val_pdf_params)

    @test typeof(p) <: Plots.Plot

    prior = get_prior(val_pdf_params)

    @test typeof(prior) <: NamedTupleDist

end

@testset "Dirichlet PDF parametrisation" begin

    local dir_pdf_params

    for i in 1:100

        K_u = rand(Uniform(2, 10))
        K_d = rand(Uniform(2, 10))
        λ_g1 = rand(Uniform(0, 1))
        λ_g2 = rand(Uniform(-1, 0))
        K_g = rand(Uniform(2, 10))
        λ_q = rand(Uniform(-1, 0))
        K_q = rand(Uniform(1, 5))
        weights = ones(9)
        θ = rand(Dirichlet(weights))

        dir_pdf_params = DirichletPDFParams(K_u=K_u, K_d=K_d, λ_g1=λ_g2,
            λ_g2=λ_g2, K_g=K_g, λ_q=λ_q, K_q=K_q,
            θ=θ)

        @test int_xtotx(dir_pdf_params) ≈ 1.0

    end

    @test typeof(dir_pdf_params) == DirichletPDFParams{Float64,Vector{Float64}}

    p = plot_input_pdfs(dir_pdf_params)

    @test typeof(p) <: Plots.Plot

    prior = get_prior(dir_pdf_params)

    @test typeof(prior) <: NamedTupleDist

end

@testset "Bernstein PDF parametrisation" begin

    local bern_pdf_params

    for i in 1:100

        U_weights = rand(Uniform(0, 100000), 4)
        D_weights = rand(Uniform(0, 100000), 4)
        λ_g1 = rand(Uniform(0, 1))
        λ_g2 = rand(Uniform(-1, 0))
        K_g = rand(Uniform(2, 10))
        λ_q = rand(Uniform(-1, 0))
        K_q = rand(Uniform(1, 5))
        weights = ones(7)
        seed = i

        bern_pdf_params = BernsteinPDFParams(U_weights=U_weights,
            D_weights=D_weights, λ_g1=λ_g2,
            λ_g2=λ_g2, K_g=K_g, λ_q=λ_q, K_q=K_q, seed=seed, weights=weights)

        @test int_xtotx(bern_pdf_params) ≈ 1.0

    end

    @test typeof(bern_pdf_params) == BernsteinPDFParams

    @test bern_pdf_params.param_type == BERNSTEIN_TYPE

    p = plot_input_pdfs(bern_pdf_params)

    @test typeof(p) <: Plots.Plot

    prior = get_prior(bern_pdf_params)

    @test typeof(prior) <: NamedTupleDist

end