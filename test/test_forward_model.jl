using PartonDensity
using Test
using Distributions, Random
using CSV
using QCDNUM

@testset "Forward model" begin

    seed = 0

    # Define different parametriations for testing
    bern_pdf_params = BernsteinPDFParams(U_weights=ones(4), D_weights=ones(4),
        λ_g1=0.4, λ_g2=-0.6,
        K_g=4.2, λ_q=-0.2, K_q=5.0,
        weights=[5.0, 5.0, 1.0, 1.0, 1.0, 0.5, 0.5])

    θ_val = [0.3027113383772114, 0.21562307667168612, 0.037706327593135316, 0.014630127386380845,
        0.03372023967525289, 0.03878412120839827, 0.007449769087935228] # Already scaled according to valence
    val_pdf_params = ValencePDFParams(λ_u=0.6, K_u=3.4,
        λ_d=0.7, K_d=4.7,
        λ_g1=0.4, λ_g2=-0.6,
        K_g=4.2, λ_q=-0.2, K_q=5.0,
        θ=θ_val)

    θ_dir = [0.23532642809197413, 0.06635605959465533, 0.2412750016688464,
        0.3719253531266372, 0.05301033816739624, 0.023468633864992115,
        0.0018706294432919719, 0.0034686372652611964, 0.0032989187769455093]
    dir_pdf_params = DirichletPDFParams(K_u=3.4, K_d=4.7,
        λ_g1=0.4, λ_g2=-0.6,
        K_g=4.2, λ_q=-0.2, K_q=5.0,
        θ=θ_dir)

    pdf_params_list = [val_pdf_params, dir_pdf_params, bern_pdf_params]

    # Initialise
    reference_data_file = string(@__DIR__, "/reference_test_data/counts_pred.csv")
    counts_pred = CSV.read(reference_data_file, NamedTuple)

    qcdnum_grid = QCDNUM.GridParams(x_min=[1.0e-3], x_weights=[1], nx=100,
        qq_bounds=[1.0e2, 3.0e4], qq_weights=[1.0, 1.0],
        nq=50, spline_interp=3)

    qcdnum_params = QCDNUM.EvolutionParams(order=2, α_S=0.118, q0=100.0,
        grid_params=qcdnum_grid, n_fixed_flav=5,
        iqc=1, iqb=1, iqt=1, weight_type=1)

    splint_params = QCDNUM.SPLINTParams(nuser=1000)
    quark_coeffs = QuarkCoefficients()

    forward_model_init(qcdnum_params, splint_params)

    # Run forward model
    for pdf_params in pdf_params_list

        counts_pred_ep, counts_pred_em = forward_model(pdf_params, qcdnum_params,
            splint_params, quark_coeffs)

        if typeof(pdf_params) == BernsteinPDFParams

            @test all(counts_pred_ep .>= 0.0)
            @test all(counts_pred_ep .<= 2.0e3)

            @test all(counts_pred_em .>= 0.0)
            @test all(counts_pred_em .<= 2.0e3)

        elseif typeof(pdf_params) == DirichletPDFParams{Float64,Vector{Float64}}

            @test all(counts_pred_ep .>= 0.0)
            @test all(counts_pred_ep .<= 1.0e3)

            @test all(counts_pred_em .>= 0.0)
            @test all(counts_pred_em .<= 1.0e3)

            @test all(counts_pred_ep .≈ counts_pred.counts_pred_ep_dir)
            @test all(counts_pred_em .≈ counts_pred.counts_pred_em_dir)

        elseif typeof(pdf_params) == ValencePDFParams{Float64,Vector{Float64}}

            @test all(counts_pred_ep .>= 0.0)
            @test all(counts_pred_ep .<= 1.0e3)

            @test all(counts_pred_em .>= 0.0)
            @test all(counts_pred_em .<= 1.0e3)

            @test all(counts_pred_ep .≈ counts_pred.counts_pred_ep_val)
            @test all(counts_pred_em .≈ counts_pred.counts_pred_em_val)

        end

        nbins = size(counts_pred_ep)[1]
        counts_obs_ep = zeros(UInt64, nbins)
        counts_obs_em = zeros(UInt64, nbins)

        for i in 1:nbins
            counts_obs_ep[i] = rand(Poisson(counts_pred_ep[i]))
            counts_obs_em[i] = rand(Poisson(counts_pred_em[i]))
        end

        sim_data = Dict{String,Any}()
        sim_data["nbins"] = nbins
        sim_data["counts_obs_ep"] = counts_obs_ep
        sim_data["counts_obs_em"] = counts_obs_em

        mktempdir() do tmp_dir

            output_file = joinpath(tmp_dir, "test_sim.h5")

            pd_write_sim(output_file, pdf_params, sim_data, MD_ZEUS_I1787035)

            new_pdf_params, new_sim_data = pd_read_sim(output_file)

            @test typeof(new_pdf_params) == typeof(pdf_params)

            @test new_sim_data == new_sim_data

        end

    end

end
