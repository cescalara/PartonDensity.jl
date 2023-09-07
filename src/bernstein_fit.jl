using ValueShapes
using ParallelProcessingTools
using DensityInterface

export get_prior, get_likelihood


function get_prior end


function get_prior(pdf_params::BernsteinPDFParams)

    prior = NamedTupleDist(
        θ_tmp=Dirichlet(pdf_params.weights),
        U_weights=[Uniform(1, 3), Uniform(1, 3), Uniform(1, 3), Uniform(1, 3)],
        D_weights=[Uniform(1, 3), Uniform(1, 3), Uniform(1, 3), Uniform(1, 3)],
        λ_g1=Uniform(0, 1),
        λ_g2=Uniform(-1, 0),
        K_g=Uniform(2, 10),
        λ_q=Uniform(-1, 0),
        K_q=Uniform(1, 5),
    )

    return prior
end

"""
Placeholder function for priors.

bspoly_params needs to be input as Vector{Int64}
since NamedTupleDist() cannot work with Vector{Vector{Int64}}
"""


function get_prior(pdf_params::BernsteinDirichletPDFParams)

    prior = NamedTupleDist(
        θ=Dirichlet(pdf_params.weights),
        initial_U=Uniform(0.0, 1.0),
        initial_D=Uniform(0.0, 1.0),
        λ_g1=Uniform(0, 1),
        λ_g2=Uniform(-1, 0),
        K_g=Uniform(2, 10),
        λ_q=Uniform(-1, 0),
        bspoly_params=[0, 4, 1, 4, 0, 5],
        K_q=Uniform(1.0, 5.0),
    )

    return prior
end


function get_likelihood(pdf_params::BernsteinPDFParams, sim_data::Dict{String,Any},
    qcdnum_params::QCDNUM.EvolutionParams, splint_params::QCDNUM.SPLINTParams,
    quark_coeffs::QuarkCoefficients, md::MetaData)

    likelihood = let d = sim_data

        counts_obs_ep = Int.(d["counts_obs_ep"])
        counts_obs_em = Int.(d["counts_obs_em"])
        nbins = d["nbins"]

        logfuncdensity(function (params)

            U_list = get_scaled_UD(Vector(params.U_weights), 2)
            D_list = get_scaled_UD(Vector(params.U_weights), 1)

            θ = get_scaled_θ(U_list, D_list, Vector(params.θ_tmp))

            pdf_params = BernsteinPDFParams(U_list=U_list, D_list=D_list,
                λ_g1=params.λ_g1, λ_g2=params.λ_g2,
                K_g=params.K_g, λ_q=params.λ_q, θ=θ, K_q=params.K_q,
                bspoly_params=Vector(params.bspoly_params)
            )

            counts_pred_ep, counts_pred_em = @critical forward_model(pdf_params, qcdnum_params,splint_params, quark_coeffs,md)

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


function get_likelihood(pdf_params::BernsteinDirichletPDFParams, sim_data::Dict{String,Any},
    qcdnum_params::QCDNUM.EvolutionParams, splint_params::QCDNUM.SPLINTParams,
    quark_coeffs::QuarkCoefficients, pos_init_u_only::Bool,md::MetaData)

    likelihood = let d = sim_data

        counts_obs_ep = Int.(d["counts_obs_ep"])
        counts_obs_em = Int.(d["counts_obs_em"])
        nbins = d["nbins"]

        logfuncdensity(function (params)

            vec_bspp = Vector(params.bspoly_params)
            bspoly_params = [[vec_bspp[Int(2 * i - 1)], vec_bspp[Int(2 * i)]] for i in 1:length(vec_bspp)/2]

            bspoly_params_d = 0

            try
                vec_bsppd = Vector(params.bspoly_params_d)
                bspoly_params_d = [[vec_bsppd[Int(2 * i - 1)], vec_bsppd[Int(2 * i)]] for i in 1:length(vec_bsppd)/2]
            catch err
                bspoly_params_d = bspoly_params
            end

            initU = Vector(params.initial_U)
            initD = Vector(params.initial_D)

            if pos_init_u_only && any(x -> x <= 0.0, initU)
                ll_value = -1000

            else

                pdf_params = BernsteinDirichletPDFParams(initial_U=initU,
                    initial_D=initD,
                    λ_g1=params.λ_g1, λ_g2=params.λ_g2,
                    K_g=params.K_g, λ_q=params.λ_q, θ=Vector(params.θ), K_q=params.K_q,
                    bspoly_params=bspoly_params,
                    bspoly_params_d=bspoly_params_d)

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

            end

            return ll_value
        end)

    end

    return likelihood
end


