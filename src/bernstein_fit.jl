using ValueShapes
using ParallelProcessingTools
using DensityInterface

export get_prior, get_likelihood
export plot_model_space, plot_data_space


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
        θ = Dirichlet(pdf_params.weights),
        initial_U = Uniform(0.,1.),
        initial_D = Uniform(0.,1.),
        λ_g1 = Uniform(0, 1),
        λ_g2 = Uniform(-1, 0),
        K_g =  Uniform(2, 10),
        λ_q = Uniform(-1, 0),
        bspoly_params = [0, 4, 1, 4, 0, 5],
        K_q=Uniform(1., 5.),
    )

    return prior
end


function get_likelihood(pdf_params::BernsteinPDFParams, sim_data::Dict{String, Any},
                        qcdnum_params::QCDNUMParameters, splint_params::SPLINTParameters,
                        quark_coeffs::QuarkCoefficients)

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
                                            bspoly_params = Vector(params.bspoly_params)
                                            )
            
            counts_pred_ep, counts_pred_em = @critical forward_model(pdf_params, qcdnum_params,
                                                                     splint_params, quark_coeffs)

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


function get_likelihood(pdf_params::BernsteinDirichletPDFParams, sim_data::Dict{String, Any},
                        qcdnum_params::QCDNUMParameters, splint_params::SPLINTParameters,
                        quark_coeffs::QuarkCoefficients)

    likelihood = let d = sim_data

        counts_obs_ep = Int.(d["counts_obs_ep"])
        counts_obs_em = Int.(d["counts_obs_em"])
        nbins = d["nbins"]

        logfuncdensity(function (params)
       
            vec_bspp = Vector(params.bspoly_params)
            bspoly_params = [[vec_bspp[Int(2*i-1)], vec_bspp[Int(2*i)]] for i in 1:length(vec_bspp)/2]

            bspoly_params_d = 0

            try
                vec_bsppd = Vector(params.bspoly_params_d)
                bspoly_params_d = [[vec_bsppd[Int(2*i-1)], vec_bsppd[Int(2*i)]] for i in 1:length(vec_bsppd)/2]
	        catch err
		        bspoly_params_d = bspoly_params
	        end
            
            initU = Vector(params.initial_U)
            initD = Vector(params.initial_D)
            
            pdf_params = BernsteinDirichletPDFParams(initial_U = initU, 
                                            initial_D = initD, 
                                            λ_g1=params.λ_g1, λ_g2=params.λ_g2,
                                            K_g=params.K_g, λ_q=params.λ_q, θ=Vector(params.θ), K_q=params.K_q,
                                            bspoly_params = bspoly_params,
                                            bspoly_params_d = bspoly_params_d)
            
            counts_pred_ep, counts_pred_em = @critical forward_model(pdf_params, qcdnum_params,
                splint_params, quark_coeffs)

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


function plot_model_space end


function plot_model_space(pdf_params::AbstractPDFParams, samples; xmin::Float64=1e-3,
    xmax::Float64=1.0, nx::Integer=50, nsamples::Integer=200,
    truth_color=:black, sample_color=:skyblue3)

    x_grid = range(xmin, stop=xmax, length=nx)

    p = plot(x_grid, [xtotx(x, pdf_params) for x in x_grid], color=truth_color,
        lw=3, label="Truth")

    sub_samples = BAT.bat_sample(samples, BAT.OrderedResampling(nsamples=nsamples)).result

    p = plot_model_space_impl(x_grid, pdf_params, sub_samples, p, color=sample_color)

    p = plot!(xaxis=:log, yaxis=:log, xlabel="x", ylabel="xtotx")
    p = ylims!(1e-5, 50.0)

    return p
end


function plot_model_space_impl(x_grid::StepRangeLen{Float64}, pdf_params::BernsteinPDFParams, samples, p; color=:skyblue3)

    for i in eachindex(samples)

        U_list = get_scaled_UD(Vector(samples.v.U_weights[i]), 2)
        D_list = get_scaled_UD(Vector(samples.v.U_weights[i]), 1)

        θ_i = get_scaled_θ(U_list, D_list, Vector(samples.v.θ_tmp[i]))
            
        vec_bspp = Vector(samples.bspoly_params[i])
        bspoly_params = [[vec_bspp[Int(2*i-1)], vec_bspp[Int(2*i)]] for i in 1:length(vec_bspp)/2]
        
        bspoly_params_d = 0

        try
            vec_bsppd = Vector(samples.bspoly_params_d[i])
            bspoly_params_d = [[vec_bsppd[Int(2*i-1)], vec_bsppd[Int(2*i)]] for i in 1:length(vec_bsppd)/2]
	    catch err
	        bspoly_params_d = bspoly_params
        end

        pdf_params_i = BernsteinPDFParams(U_list=U_list, D_list=D_list,
            λ_g1=samples.v.λ_g1[i], λ_g2=samples.v.λ_g2[i],
            K_g=samples.v.K_g[i], λ_q=samples.v.λ_q[i], K_q=samples.v.K_q[i],
            θ=θ_i,
            bspoly_params = bspoly_params,
            bspoly_params_d = bspoly_params_d)
        
        p = plot!(x_grid, [xtotx(x, pdf_params_i) for x in x_grid], color=color, lw=3,
            alpha=0.01, label="")

    end

    return p
end


function plot_model_space_impl(x_grid::StepRangeLen{Float64}, pdf_params::BernsteinDirichletPDFParams, samples, p; color=:skyblue3)

    for i in eachindex(samples)
        
        vec_bspp = Vector(samples.bspoly_params[i])
        bspoly_params = [[vec_bspp[Int(2*i-1)], vec_bspp[Int(2*i)]] for i in 1:length(vec_bspp)/2]
        
        bspoly_params_d = 0

        try
            vec_bsppd = Vector(samples.bspoly_params_d[i])
            bspoly_params_d = [[vec_bsppd[Int(2*i-1)], vec_bsppd[Int(2*i)]] for i in 1:length(vec_bsppd)/2]
	    catch err
	        bspoly_params_d = bspoly_params
        end

        pdf_params_i = BernsteinDirichletPDFParams(initial_U=[samples.v.initial_U[i]], initial_D=[samples.v.initial_D[i]],
                                        λ_g1=samples.v.λ_g1[i], λ_g2=samples.v.λ_g2[i],
                                        K_g=samples.v.K_g[i], λ_q=samples.v.λ_q[i], 
                                        θ=Vector(samples.v.θ[i]),
                                        bspoly_params = bspoly_params,
                                        bspoly_params_d = bspoly_params_d)
        
        p = plot!(x_grid, [xtotx(x, pdf_params_i) for x in x_grid], color=color, lw=3,
                  alpha=0.01, label="")
        
    end

    return p
end


function plot_data_space end


function plot_data_space(pdf_params::AbstractPDFParams, sim_data::Dict{String,Any}, samples,
    qcdnum_grid::QCDNUMGrid, qcdnum_params::QCDNUMParameters,
    splint_params::SPLINTParameters, quark_coeffs::QuarkCoefficients;
    ep_color=:firebrick, em_color=:teal, nsamples::Integer=100, plot_size=(1000, 500))

    forward_model_init(qcdnum_grid, qcdnum_params, splint_params)

    nbins = sim_data["nbins"]

    p1 = scatter(1:nbins, sim_data["counts_obs_ep"], label="Observed counts (eP)",
        color=ep_color, lw=3)
    p2 = scatter(1:nbins, sim_data["counts_obs_em"], label="Observed counts (eM)",
        color=em_color, lw=3)

    sub_samples = BAT.bat_sample(samples, BAT.OrderedResampling(nsamples=nsamples)).result

    p1, p2 = plot_data_space_impl(pdf_params, sub_samples, qcdnum_params,
        splint_params, quark_coeffs, p1, p2, nbins)

    plot(p1, p2, size=plot_size, xlabel="Bin number")
end


function plot_data_space_impl(pdf_params::BernsteinPDFParams, samples, qcdnum_params::QCDNUMParameters,
    splint_params::SPLINTParameters, quark_coeffs::QuarkCoefficients,
    p1, p2, nbins::Integer; ep_color=:firebrick, em_color=:teal)

    for i in eachindex(samples)

        counts_obs_ep_i = zeros(UInt64, nbins)
        counts_obs_em_i = zeros(UInt64, nbins)
        
        U_list = get_scaled_UD(Vector(samples.v.U_weights[i]), 2)
        D_list = get_scaled_UD(Vector(samples.v.U_weights[i]), 1)

        θ_i = get_scaled_θ(U_list, D_list, Vector(samples.v.θ_tmp[i]))

        vec_bspp = Vector(samples.bspoly_params[i])
        bspoly_params = [[vec_bspp[Int(2*i-1)], vec_bspp[Int(2*i)]] for i in 1:length(vec_bspp)/2]
        
        bspoly_params_d = 0

        try
            vec_bsppd = Vector(samples.bspoly_params_d[i])
            bspoly_params_d = [[vec_bsppd[Int(2*i-1)], vec_bsppd[Int(2*i)]] for i in 1:length(vec_bsppd)/2]
	    catch err
	        bspoly_params_d = bspoly_params
        end

        pdf_params_i = BernsteinPDFParams(U_list=U_list, D_list=D_list,
            λ_g1=samples.v.λ_g1[i], λ_g2=samples.v.λ_g2[i],
            K_g=samples.v.K_g[i], λ_q=samples.v.λ_q[i], K_q=samples.v.K_q[i],
            θ=θ_i, bspoly_params = bspoly_params,
            bspoly_params_d = bspoly_params_d)

        counts_pred_ep_i, counts_pred_em_i = forward_model(pdf_params_i, qcdnum_params,
            splint_params, quark_coeffs)

        for j in 1:nbins

            if counts_pred_ep_i[j] < 0

                @warn "Predicted counts (eP) found to be < 0, setting to 0" i j counts_pred_ep_i[j]
                counts_pred_ep_i[j] = 0

            end

            if counts_pred_em_i[j] < 0

                @warn "Predicted counts (eM) found to be < 0, setting to 0" i j counts_pred_em_i[j]
                counts_pred_em_i[j] = 0

            end

            counts_obs_ep_i[j] = rand(Poisson(counts_pred_ep_i[j]))
            counts_obs_em_i[j] = rand(Poisson(counts_pred_em_i[j]))

        end

        p1 = scatter!(p1, 1:nbins, counts_obs_ep_i, label="", color=ep_color,
            lw=3, alpha=0.01)
        p2 = scatter!(p2, 1:nbins, counts_obs_em_i, label="", color=em_color,
            lw=3, alpha=0.01)

    end

    return p1, p2
end


function plot_data_space_impl(pdf_params::BernsteinDirichletPDFParams, samples, qcdnum_params::QCDNUMParameters,
                              splint_params::SPLINTParameters, quark_coeffs::QuarkCoefficients,
                              p1, p2, nbins::Integer; ep_color=:firebrick, em_color=:teal)

    for i in eachindex(samples)
        
        counts_obs_ep_i = zeros(UInt64, nbins)
        counts_obs_em_i = zeros(UInt64, nbins)
            
        vec_bspp = Vector(samples.bspoly_params[i])
        bspoly_params = [[vec_bspp[Int(2*i-1)], vec_bspp[Int(2*i)]] for i in 1:length(vec_bspp)/2]
        
        bspoly_params_d = 0
        
        try
            vec_bsppd = Vector(samples.bspoly_params_d[i])
            bspoly_params_d = [[vec_bsppd[Int(2*i-1)], vec_bsppd[Int(2*i)]] for i in 1:length(vec_bsppd)/2]
	    catch err
	        bspoly_params_d = bspoly_params
        end

        pdf_params_i = BernsteinDirichletPDFParams(initial_U=[samples.v.initial_U[i]], initial_D=[samples.v.initial_D[i]],
                                        λ_g1=samples.v.λ_g1[i], λ_g2=samples.v.λ_g2[i],
                                        K_g=samples.v.K_g[i], λ_q=samples.v.λ_q[i], 
                                        θ=Vector(samples.v.θ[i]),
                                        bspoly_params = bspoly_params,
                                        bspoly_params_d = bspoly_params_d)
        
        counts_pred_ep_i, counts_pred_em_i = forward_model(pdf_params_i, qcdnum_params, 
                                                           splint_params, quark_coeffs)
        
        for j in 1:nbins

            if counts_pred_ep_i[j] < 0

                @warn "Predicted counts (eP) found to be < 0, setting to 0" i j counts_pred_ep_i[j]
                counts_pred_ep_i[j] = 0
                
            end

            if counts_pred_em_i[j] < 0

                @warn "Predicted counts (eM) found to be < 0, setting to 0" i j counts_pred_em_i[j]
                counts_pred_em_i[j] = 0
        
            end
            
            counts_obs_ep_i[j] = rand(Poisson(counts_pred_ep_i[j]))
            counts_obs_em_i[j] = rand(Poisson(counts_pred_em_i[j]))
            
        end
        
        p1 = scatter!(p1, 1:nbins, counts_obs_ep_i, label="", color=ep_color, 
                      lw=3, alpha=0.01)
        p2 = scatter!(p2, 1:nbins, counts_obs_em_i, label="", color=em_color, 
                      lw=3, alpha=0.01)
    
    end
  
    return p1, p2
end
