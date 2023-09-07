using ValueShapes
using ParallelProcessingTools
using DensityInterface

export get_prior, get_likelihood

"""
    get_prior(pdf_params)

Get a suitable prior for a given PDF parametrisation.
"""
function get_prior end


function get_prior(pdf_params::ValencePDFParams)

    prior = NamedTupleDist(
        θ_tmp=Dirichlet(ones(7)),
        λ_u=Uniform(0, 1),
        K_u=Uniform(2, 10),
        λ_d=Uniform(0, 1),
        K_d=Uniform(2, 10),
        λ_g1=Uniform(0, 1),
        λ_g2=Uniform(-1, 0),
        K_g=Uniform(2, 10),
        λ_q=Uniform(-1, 0),
        K_q=Uniform(1, 5),
    )

    return prior
end


function get_prior(pdf_params::DirichletPDFParams)

    prior = NamedTupleDist(
        θ=Dirichlet(ones(9)),
        K_u=Uniform(2, 10),
        K_d=Uniform(2, 10),
        λ_g1=Uniform(0, 1),
        λ_g2=Uniform(-1, 0),
        K_g=Uniform(2, 10),
        λ_q=Uniform(-1, 0),
        K_q=Uniform(1, 5),
    )

    return prior
end


"""
    get_likelihood(pdf_params)

Get a suitable likelihood for a given PDF parametrisation.
"""
function get_likelihood end


function get_likelihood(pdf_params::ValencePDFParams, sim_data::Dict{String,Any},
    qcdnum_params::QCDNUM.EvolutionParams, splint_params::QCDNUM.SPLINTParams,
    quark_coeffs::QuarkCoefficients)

    likelihood = let d = sim_data

        counts_obs_ep = Int.(d["counts_obs_ep"])
        counts_obs_em = Int.(d["counts_obs_em"])
        nbins = d["nbins"]

        logfuncdensity(function (params)

            θ = get_scaled_θ(params.λ_u, params.K_u, params.λ_d,
                params.K_d, Vector(params.θ_tmp))

            pdf_params = ValencePDFParams(λ_u=params.λ_u, K_u=params.K_u, λ_d=params.λ_d,
                K_d=params.K_d, λ_g1=params.λ_g1, λ_g2=params.λ_g2,
                K_g=params.K_g, λ_q=params.λ_q, K_q=params.K_q, θ=θ)

            counts_pred_ep, counts_pred_em = @critical forward_model(pdf_params, qcdnum_params,splint_params, quark_coeffs,MD_ZEUS_I1787035)

            ll_value = 0.0
            for i in 1:nbins

                if counts_pred_ep[i] < 0
                    @debug "counts_pred_ep[i] < 0, setting to 0" i counts_pred_ep[i]
                    counts_pred_ep[i] = 0
                end

                if counts_pred_em[i] < 0
                    @debug "counts_pred_em[i] < 0, setting to 0" i counts_pred_em[i]
                    counts_pred_em[i] = 0
                end

                ll_value += logpdf(Poisson(counts_pred_ep[i]), counts_obs_ep[i])
                ll_value += logpdf(Poisson(counts_pred_em[i]), counts_obs_em[i])
            end

            return ll_value
        end)

    end

    return likelihood
end


function get_likelihood(pdf_params::DirichletPDFParams, sim_data::Dict{String,Any},
    qcdnum_params::QCDNUM.EvolutionParams, splint_params::QCDNUM.SPLINTParams,
    quark_coeffs::QuarkCoefficients)

    likelihood = let d = sim_data

        counts_obs_ep = Int.(d["counts_obs_ep"])
        counts_obs_em = Int.(d["counts_obs_em"])
        nbins = d["nbins"]

        logfuncdensity(function (params)

            pdf_params = DirichletPDFParams(K_u=params.K_u, K_d=params.K_d,
                λ_g1=params.λ_g1, λ_g2=params.λ_g2,
                K_g=params.K_g, λ_q=params.λ_q, K_q=params.K_q,
                θ=Vector(params.θ))

            # Ensure u-valence weight > d-valence weight
            if params.θ[2] > params.θ[1]

                return -Inf

            end

            counts_pred_ep, counts_pred_em = @critical forward_model(pdf_params, qcdnum_params,splint_params, quark_coeffs,MD_ZEUS_I1787035)

            ll_value = 0.0
            for i in 1:nbins

                if counts_pred_ep[i] < 0
                    @debug "counts_pred_ep[i] < 0, setting to 0" i counts_pred_ep[i]
                    counts_pred_ep[i] = 0
                end

                if counts_pred_em[i] < 0
                    @debug "counts_pred_em[i] < 0, setting to 0" i counts_pred_em[i]
                    counts_pred_em[i] = 0
                end

                ll_value += logpdf(Poisson(counts_pred_ep[i]), counts_obs_ep[i])
                ll_value += logpdf(Poisson(counts_pred_em[i]), counts_obs_em[i])
            end

            return ll_value
        end)

    end

    return likelihood
end
