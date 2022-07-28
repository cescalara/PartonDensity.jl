using PartonDensity
using Test
using Random, Distributions
using Plots, ValueShapes
import BAT

@testset "Bernstein PDF parametrisation" begin

    local bern_pdf_params
    
    for i in 1:100

        U_weights = rand(Uniform(0, 100000), 4)
        D_weights = rand(Uniform(0, 100000), 4)
        λ_g1 = rand(Uniform(0, 1))
        λ_g2 = rand(Uniform(-1, 0))
        K_g = rand(Uniform(2, 10))
        λ_q = rand(Uniform(-1, 0))
        weights = ones(7)
        
        bern_pdf_params = BernsteinPDFParams(U_weights = U_weights,
                                          D_weights = D_weights, λ_g1=λ_g2,
                                          λ_g2=λ_g2, K_g=K_g, λ_q=λ_q,
                                          weights=weights)

        
        @test int_xtotx(bern_pdf_params) ≈ 1.0
          
    end

    @test typeof(bern_pdf_params) == BernsteinPDFParams

    @test bern_pdf_params.param_type == BERNSTEIN_TYPE

    p = plot_input_pdfs(bern_pdf_params)
    
    @test typeof(p) <: Plots.Plot

    prior = get_prior(bern_pdf_params)
    
    @test typeof(prior) <: NamedTupleDist 
    
end


