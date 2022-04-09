using PartonDensity
using Test
using Random, Distributions
using Plots

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
        weights = ones(7)
        
        val_pdf_params = ValencePDFParams(λ_u=λ_u, K_u=K_u,
                                          λ_d=λ_d, K_d=K_d, λ_g1=λ_g2,
                                          λ_g2=λ_g2, K_g=K_g, λ_q=λ_q,
                                          weights=weights)

        
        @test int_xtotx(val_pdf_params) ≈ 1.0
          
    end

    @test typeof(val_pdf_params) == ValencePDFParams

    @test val_pdf_params.param_type == VALENCE_TYPE

    p = plot_input_pdfs(val_pdf_params)
    
    @test typeof(p) <: Plots.Plot
    
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
        weights = ones(9)
        
        dir_pdf_params = DirichletPDFParams(K_u=K_u, K_d=K_d, λ_g1=λ_g2,
                                            λ_g2=λ_g2, K_g=K_g, λ_q=λ_q,
                                            weights=weights)

        @test int_xtotx(dir_pdf_params) ≈ 1.0   
          
    end

    @test typeof(dir_pdf_params) == DirichletPDFParams

    @test dir_pdf_params.param_type == DIRICHLET_TYPE

    p = plot_input_pdfs(dir_pdf_params)

    @test typeof(p) <: Plots.Plot
    
end
