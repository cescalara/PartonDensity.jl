using BAT
using ValueShapes
using ParallelProcessingTools
using DensityInterface

export get_prior, get_likelihood
export plot_model_space, plot_data_space

"""
    get_prior(pdf_params)

Get a suitable prior for a given PDF parametrisation.
"""
function get_prior end


function get_prior(pdf_params::ValencePDFParams)

    prior = NamedTupleDist(
        θ_tmp = Dirichlet(pdf_params.weights),
        λ_u = Uniform(0, 1), 
        K_u = Uniform(2, 10),
        λ_d = Uniform(0, 1), 
        K_d = Uniform(2, 10),
        λ_g1 = Uniform(0, 1),
        λ_g2 = Uniform(-1, 0),
        K_g =  Uniform(2, 10),
        λ_q = Uniform(-1, 0),
    )

    return prior
end


function get_prior(pdf_params::DirichletPDFParams)

    prior = NamedTupleDist(
        θ = Dirichlet(pdf_params.weights),
        K_u = Uniform(2, 10),
        K_d = Uniform(2, 10),
        λ_g1 = Uniform(0, 1),
        λ_g2 = Uniform(-1, 0),
        K_g =  Uniform(2, 10),
        λ_q = Uniform(-1, 0),
    )

    return prior
end


"""
    get_likelihood(pdf_params)

Get a suitable likelihood for a given PDF parametrisation.
"""
function get_likelihood end


function get_likelihood(pdf_params::ValencePDFParams, sim_data::Dict{String, Any},
                        qcdnum_params::QCDNUMParameters, splint_params::SPLINTParameters,
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
                                          K_g=params.K_g, λ_q=params.λ_q, θ=θ)
            
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


function get_likelihood(pdf_params::DirichletPDFParams, sim_data::Dict{String, Any},
                        qcdnum_params::QCDNUMParameters, splint_params::SPLINTParameters,
                        quark_coeffs::QuarkCoefficients)

    likelihood = let d = sim_data

        counts_obs_ep = Int.(d["counts_obs_ep"])
        counts_obs_em = Int.(d["counts_obs_em"])
        nbins = d["nbins"]

        logfuncdensity(function (params)
            
            pdf_params = DirichletPDFParams(K_u=params.K_u, K_d=params.K_d,
                                            λ_g1=params.λ_g1, λ_g2=params.λ_g2,
                                            K_g=params.K_g, λ_q=params.λ_q,
                                            θ=Vector(params.θ))

            # Ensure u-valence weight > d-valence weight
            if params.θ[2] > params.θ[1]
                   
                return -Inf
            
            end
                     
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


"""
    plot_model_space(pdf_params, samples)

Compare truth and posterior samples in the model space.
"""
function plot_model_space end


function plot_model_space(pdf_params::AbstractPDFParams, samples; xmin::Float64=1e-3,
                          xmax::Float64=1.0, nx::Integer=50, nsamples::Integer=200,
                          truth_color=:black, sample_color=:skyblue3)

    x_grid = range(xmin, stop=xmax, length=nx)

    p = plot(x_grid, [xtotx(x, pdf_params) for x in x_grid], color=truth_color,
             lw=3, label="Truth")

    sub_samples = bat_sample(samples, OrderedResampling(nsamples=nsamples)).result

    p = plot_model_space_impl(x_grid, pdf_params, sub_samples, p, color=sample_color)

    p = plot!(xaxis=:log, yaxis=:log, xlabel="x", ylabel="xtotx")
    p = ylims!(1e-5, 50.0)

    return p
end

    
function plot_model_space_impl(x_grid::StepRangeLen{Float64}, pdf_params::ValencePDFParams, samples, p; color=:skyblue3)

    for i in eachindex(samples)

        θ_i = get_scaled_θ(samples.v.λ_u[i], samples.v.K_u[i], samples.v.λ_d[i],
                           samples.v.K_d[i], Vector(samples.v.θ_tmp[i]))
        
        pdf_params_i = ValencePDFParams(λ_u=samples.v.λ_u[i], K_u=samples.v.K_u[i], 
                                        λ_d=samples.v.λ_d[i], K_d=samples.v.K_d[i],
                                        λ_g1=samples.v.λ_g1[i], λ_g2=samples.v.λ_g2[i],
                                        K_g=samples.v.K_g[i], λ_q=samples.v.λ_q[i], 
                                        θ=θ_i)
        p = plot!(x_grid, [xtotx(x, pdf_params_i) for x in x_grid], color=color, lw=3,
                  alpha=0.01, label="")
        
    end

    return p
end


function plot_model_space_impl(x_grid::StepRangeLen{Float64}, pdf_params::DirichletPDFParams, samples, p; color=:skyblue3)

    for i in eachindex(samples)

        pdf_params_i = DirichletPDFParams(K_u=samples.v.K_u[i], K_d=samples.v.K_d[i],
                                          λ_g1=samples.v.λ_g1[i], λ_g2=samples.v.λ_g2[i],
                                          K_g=samples.v.K_g[i], λ_q=samples.v.λ_q[i], 
                                          θ=Vector(samples.v.θ[i]))
        p = plot!(x_grid, [xtotx(x, pdf_params_i) for x in x_grid], color=color, lw=3,
                  alpha=0.01, label="")
        
    end

    return p
end


"""
    plot_data_space(pdf_params, sim_data, samples, qcdnum_grid, 
                    qcdnum_params, splint_params, quark_coeffs)

Compare truth and posterior samples in the data space.
"""
function plot_data_space end


function plot_data_space(pdf_params::AbstractPDFParams, sim_data::Dict{String, Any}, samples,
                         qcdnum_grid::QCDNUMGrid, qcdnum_params::QCDNUMParameters,
                         splint_params::SPLINTParameters, quark_coeffs::QuarkCoefficients;
                         ep_color=:firebrick, em_color=:teal, nsamples::Integer=100, plot_size=(1000, 500))

    forward_model_init(qcdnum_grid, qcdnum_params, splint_params)

    nbins = sim_data["nbins"]

    p1 = scatter(1:nbins, sim_data["counts_obs_ep"], label="Observed counts (eP)", 
                 color=ep_color, lw=3)
    p2 = scatter(1:nbins, sim_data["counts_obs_em"], label="Observed counts (eM)", 
                 color=em_color, lw=3)

    sub_samples = bat_sample(samples, OrderedResampling(nsamples=nsamples)).result

    p1, p2 = plot_data_space_impl(pdf_params, sub_samples, qcdnum_params,
                                  splint_params, quark_coeffs, p1, p2, nbins)
      
    plot(p1, p2, size=plot_size, xlabel="Bin number")
end


function plot_data_space_impl(pdf_params::ValencePDFParams, samples, qcdnum_params::QCDNUMParameters,
                              splint_params::SPLINTParameters, quark_coeffs::QuarkCoefficients,
                              p1, p2, nbins::Integer; ep_color=:firebrick, em_color=:teal)

    for i in eachindex(samples)
        
        counts_obs_ep_i = zeros(UInt64, nbins)
        counts_obs_em_i = zeros(UInt64, nbins)

        θ_i = get_scaled_θ(samples.v.λ_u[i], samples.v.K_u[i], samples.v.λ_d[i],
                           samples.v.K_d[i], Vector(samples.v.θ_tmp[i]))


        pdf_params_i = ValencePDFParams(λ_u=samples.v.λ_u[i], K_u=samples.v.K_u[i], 
                                        λ_d=samples.v.λ_d[i], K_d=samples.v.K_d[i],
                                        λ_g1=samples.v.λ_g1[i], λ_g2=samples.v.λ_g2[i],
                                        K_g=samples.v.K_g[i], λ_q=samples.v.λ_q[i], 
                                        θ=θ_i)
        
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


function plot_data_space_impl(pdf_params::DirichletPDFParams, samples, qcdnum_params::QCDNUMParameters,
                              splint_params::SPLINTParameters, quark_coeffs::QuarkCoefficients,
                              p1, p2, nbins::Integer; ep_color=:firebrick, em_color=:teal)

    for i in eachindex(samples)
        
        counts_obs_ep_i = zeros(UInt64, nbins)
        counts_obs_em_i = zeros(UInt64, nbins)

        pdf_params_i = DirichletPDFParams(K_u=samples.v.K_u[i], K_d=samples.v.K_d[i],
                                          λ_g1=samples.v.λ_g1[i], λ_g2=samples.v.λ_g2[i],
                                          K_g=samples.v.K_g[i], λ_q=samples.v.λ_q[i], 
                                          θ=Vector(samples.v.θ[i]))
        
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

