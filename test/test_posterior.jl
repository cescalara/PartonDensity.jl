using PartonDensity
using BAT, DensityInterface, ParallelProcessingTools
using Random, Distributions, ValueShapes
using Test

@testset "DirichletPDFParams simulation + posterior check" begin

    rng = MersenneTwister(42)

    # Simulate from the DirichletPDFParams model
    weights = [30.0, 15.0, 12.0, 6.0, 3.6, 0.85, 0.85, 0.85, 0.85]
    θ = rand(rng, Dirichlet(weights))
    pdf_params = DirichletPDFParams(K_u=4.0, K_d=4.0, λ_g1=1.5, λ_g2=-0.4, K_g=6.0,
        λ_q=-0.25, K_q=5.0, θ=θ)

    # first specify QCDNUM inputs
    qcdnum_grid = QCDNUM.GridParams(x_min=[1.0e-3, 1.0e-1, 5.0e-1], x_weights=[1, 2, 2], nx=100,
        qq_bounds=[1.0e2, 3.0e4], qq_weights=[1.0, 1.0], nq=50, spline_interp=3)
    qcdnum_params = QCDNUM.EvolutionParams(order=2, α_S=0.118, q0=100.0, grid=qcdnum_grid,
        n_fixed_flav=5, iqc=1, iqb=1, iqt=1, weight_type=1)

    # now SPLINT and quark coefficients
    splint_params = QCDNUM.SPLINTParams()
    quark_coeffs = QuarkCoefficients()

    # initialise QCDNUM
    forward_model_init(qcdnum_params, splint_params)

    # run forward model 
    counts_pred_ep, counts_pred_em = forward_model(pdf_params, qcdnum_params,
        splint_params, quark_coeffs)

    # take a poisson sample
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

    # Define posterior
    prior = NamedTupleDist(
        θ=Dirichlet(weights),
        K_u=Truncated(Normal(3.5, 0.5), 2.0, 5.0),
        K_d=Truncated(Normal(3.5, 0.5), 2.0, 5.0),
        λ_g1=Uniform(1.0, 2.0),
        λ_g2=Uniform(-0.5, -0.1),
        K_g=Truncated(Normal(4.0, 1.5), 2.0, 5.0),
        λ_q=Uniform(-0.5, -0.1),
        K_q=Truncated(Normal(5.0, 1.5), 3.0, 10.0),
    )

    likelihood = let d = sim_data

        counts_obs_ep = d["counts_obs_ep"]
        counts_obs_em = d["counts_obs_em"]
        nbins = d["nbins"]

        logfuncdensity(function (params)

            pdf_params = DirichletPDFParams(K_u=params.K_u, K_d=params.K_d, λ_g1=params.λ_g1, λ_g2=params.λ_g2,
                K_g=params.K_g, λ_q=params.λ_q, K_q=params.K_q, θ=Vector(params.θ))

            #Ensure u-valence weight > d-valence weight
            if params.θ[2] > params.θ[1]

                return -Inf

            end

            counts_pred_ep, counts_pred_em = @critical forward_model(pdf_params,
                qcdnum_params, splint_params, quark_coeffs)

            ll_value = 0.0
            for i in 1:nbins
                ll_value += logpdf(Poisson(counts_pred_ep[i]), counts_obs_ep[i])
                ll_value += logpdf(Poisson(counts_pred_em[i]), counts_obs_em[i])
            end

            return ll_value
        end)
    end

    posterior = PosteriorDensity(likelihood, prior)

    # Evaluate posterior for random samples of the prior
    for i in 1:100

        rng = MersenneTwister(i)
        log_density = BAT.checked_logdensityof(posterior, rand(rng, prior))

        @test typeof(log_density) == Float64
        @test !isnan(log_density)

    end

end

@testset "ValencePDFParams simulation + posterior check" begin

    rng = MersenneTwister(42)

    # Simulate from the ValencePDFParams model
    weights = [5.0, 5.0, 1.0, 1.0, 1.0, 0.5, 0.5]
    λ_u = 0.64
    K_u = 3.38
    λ_d = 0.67
    K_d = 4.73

    θ_tmp = rand(rng, Dirichlet(weights))
    θ = get_scaled_θ(λ_u, K_u, λ_d, K_d, θ_tmp)
    pdf_params = ValencePDFParams(λ_u=λ_u, K_u=K_u, λ_d=λ_d, K_d=K_d,
        λ_g1=0.50, λ_g2=-0.63, K_g=4.23, λ_q=-0.23, K_q=5.0, θ=θ)

    # first specify QCDNUM inputs
    qcdnum_grid = QCDNUM.GridParams(x_min=[1.0e-3, 1.0e-1, 5.0e-1], x_weights=[1, 2, 2], nx=100,
        qq_bounds=[1.0e2, 3.0e4], qq_weights=[1.0, 1.0], nq=50, spline_interp=3)
    qcdnum_params = QCDNUM.EvolutionParams(order=2, α_S=0.118, q0=100.0, grid=qcdnum_grid,
        n_fixed_flav=5, iqc=1, iqb=1, iqt=1, weight_type=1)

    # now SPLINT and quark coefficients
    splint_params = QCDNUM.SPLINTParams()
    quark_coeffs = QuarkCoefficients()

    # initialise QCDNUM
    forward_model_init(qcdnum_params, splint_params)

    # run forward model 
    counts_pred_ep, counts_pred_em = forward_model(pdf_params, qcdnum_params,
        splint_params, quark_coeffs)

    # take a poisson sample
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

    # Define posterior
    prior = NamedTupleDist(
        θ_tmp=Dirichlet(weights),
        λ_u=Truncated(Normal(pdf_params.λ_u, 1), 0, 1),
        K_u=Truncated(Normal(3.5, 0.5), 2.0, 5.0),
        λ_d=Truncated(Normal(pdf_params.λ_d, 1), 0, 1),
        K_d=Truncated(Normal(3.5, 0.5), 2.0, 5.0),
        λ_g1=Uniform(1.0, 2.0),
        λ_g2=Uniform(-0.5, -0.1),
        K_g=Truncated(Normal(4.0, 1.5), 2.0, 5.0),
        λ_q=Uniform(-0.5, -0.1),
        K_q=Truncated(Normal(5.0, 1.5), 3.0, 10.0),
    )

    likelihood = let d = sim_data

        counts_obs_ep = d["counts_obs_ep"]
        counts_obs_em = d["counts_obs_em"]
        nbins = d["nbins"]

        logfuncdensity(function (params)

            θ = get_scaled_θ(params.λ_u, params.K_u, params.λ_d, params.K_d, Vector(params.θ_tmp))
            pdf_params = ValencePDFParams(λ_u=params.λ_u, K_u=params.K_u, λ_d=params.λ_d, K_d=params.K_d, λ_g1=params.λ_g1, λ_g2=params.λ_g2,
                K_g=params.K_g, λ_q=params.λ_q, K_q=params.K_q, θ=θ)

            counts_pred_ep, counts_pred_em = @critical forward_model(pdf_params,
                qcdnum_params, splint_params, quark_coeffs)

            ll_value = 0.0
            for i in 1:nbins
                ll_value += logpdf(Poisson(counts_pred_ep[i]), counts_obs_ep[i])
                ll_value += logpdf(Poisson(counts_pred_em[i]), counts_obs_em[i])
            end

            return ll_value
        end)
    end

    posterior = PosteriorDensity(likelihood, prior)

    # Evaluate posterior for random samples of the prior
    for i in 1:100

        rng = MersenneTwister(i)
        log_density = BAT.checked_logdensityof(posterior, rand(rng, prior))

        @test typeof(log_density) == Float64
        @test !isnan(log_density)

    end

end