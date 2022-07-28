using PartonDensity
using Test
using Distributions, Random

@testset "Forward model" begin

    # Define different parametriations for testing
    bern_pdf_params = BernsteinPDFarams(U_weights = ones(4), D_weights = ones(4),
                                      λ_g1=0.4, λ_g2=-0.6,
                                      K_g=4.2, λ_q=-0.2, 
                                      weights=[5., 5., 1., 1., 1., 0.5, 0.5])
    
    val_pdf_params = ValencePDFParams(λ_u=0.6, K_u=3.4,
                                      λ_d=0.7, K_d=4.7,
                                      λ_g1=0.4, λ_g2=-0.6,
                                      K_g=4.2, λ_q=-0.2, 
                                      weights=[5., 5., 1., 1., 1., 0.5, 0.5])

    dir_pdf_params = DirichletPDFParams(K_u=3.4, K_d=4.7,
                                        λ_g1=0.4, λ_g2=-0.6,
                                        K_g=4.2, λ_q=-0.2,
                                        weights=[3., 1., 5., 5., 1., 1., 1., 0.5, 0.5])

    pdf_params_list = [val_pdf_params, dir_pdf_params, bern_pdf_params]


    # Initialise
    qcdnum_grid = QCDNUMGrid(x_min=[1.0e-3], x_weights=[1], nx=100,
                             qq_bounds=[1.0e2, 3.0e4], qq_weights=[1.0, 1.0],
                             nq=50, spline_interp=3)
    
    qcdnum_params = QCDNUMParameters(order=2, α_S=0.118, q0=100.0,
                                     grid=qcdnum_grid, n_fixed_flav=5,
                                     iqc=1, iqb=1, iqt=1, weight_type=1)

    splint_params = SPLINTParameters(nuser=1000)
    quark_coeffs = QuarkCoefficients()

    forward_model_init(qcdnum_grid, qcdnum_params, splint_params)

    # Run forward model
    for pdf_params in pdf_params_list

        counts_pred_ep, counts_pred_em = forward_model(pdf_params, qcdnum_params, 
                                                       splint_params, quark_coeffs)

        @test all(counts_pred_ep .>= 0.0)
        @test all(counts_pred_ep .<= 1.0e3)

        @test all(counts_pred_em .> 0.0)
        @test all(counts_pred_em .<= 1.0e3)

        nbins = size(counts_pred_ep)[1]
        counts_obs_ep = zeros(UInt64, nbins)
        counts_obs_em = zeros(UInt64, nbins)

        for i in 1:nbins
            counts_obs_ep[i] = rand(Poisson(counts_pred_ep[i]))
            counts_obs_em[i] = rand(Poisson(counts_pred_em[i]))
        end

        sim_data = Dict{String, Any}()
        sim_data["nbins"] = nbins
        sim_data["counts_obs_ep"] = counts_obs_ep
        sim_data["counts_obs_em"] = counts_obs_em

        mktempdir() do tmp_dir
        
            output_file = joinpath(tmp_dir, "test_sim.h5")

            pd_write_sim(output_file, pdf_params, sim_data)

            new_pdf_params, new_sim_data = pd_read_sim(output_file)

            @test typeof(new_pdf_params) == typeof(pdf_params)
            
            @test new_sim_data == new_sim_data
        
        end
   
    end

end
