using BAT

export plot_model_space, plot_data_space

"""
    plot_model_space(pdf_params, samples)

Compare truth and posterior samples in the model space.
"""
function plot_model_space(pdf_params::PDFParameters, samples; xmin::Float64=1e-3,
                          xmax::Float64=1.0, nx::Integer=50, nsamples::Integer=200,
                          truth_color=:black, sample_color=:skyblue3)

    x_grid = range(xmin, stop=xmax, length=nx)

    p = plot(x_grid, [xtotx(x, pdf_params) for x in x_grid], color=truth_color,
             lw=3, label="Truth")

    sub_samples = bat_sample(samples, OrderedResampling(nsamples=nsamples)).result
    
    for i in eachindex(sub_samples)

        θ_i = get_scaled_θ(sub_samples.v.λ_u[i], sub_samples.v.K_u[i], sub_samples.v.λ_d[i],
                           sub_samples.v.K_d[i], Vector(sub_samples.v.θ_tmp[i]))
        
        pdf_params_i = PDFParameters(λ_u=sub_samples.v.λ_u[i], K_u=sub_samples.v.K_u[i], 
                                     λ_d=sub_samples.v.λ_d[i], K_d=sub_samples.v.K_d[i],
                                     λ_g1=sub_samples.v.λ_g1[i], λ_g2=sub_samples.v.λ_g2[i],
                                     K_g=sub_samples.v.K_g[i], λ_q=sub_samples.v.λ_q[i], 
                                     θ=θ_i)
        p = plot!(x_grid, [xtotx(x, pdf_params_i) for x in x_grid], color=sample_color, lw=3,
                  alpha=0.01, label="")
        
    end

    p = plot!(xaxis=:log, yaxis=:log, xlabel="x", ylabel="xtotx")
    p = ylims!(1e-5, 50.0)

    return p
end

"""
    plot_data_space(pdf_params, sim_data, samples, qcdnum_grid, 
                    qcdnum_params, splint_params, quark_coeffs)

Compare truth and posterior samples in the data space.
"""
function plot_data_space(pdf_params::PDFParameters, sim_data::Dict{String, Any}, samples,
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
    
    for i in eachindex(sub_samples)
        
        counts_obs_ep_i = zeros(UInt64, nbins)
        counts_obs_em_i = zeros(UInt64, nbins)

        θ_i = get_scaled_θ(sub_samples.v.λ_u[i], sub_samples.v.K_u[i], sub_samples.v.λ_d[i],
                           sub_samples.v.K_d[i], Vector(sub_samples.v.θ_tmp[i]))


        pdf_params_i = PDFParameters(λ_u=sub_samples.v.λ_u[i], K_u=sub_samples.v.K_u[i], 
                                     λ_d=sub_samples.v.λ_d[i], K_d=sub_samples.v.K_d[i],
                                     λ_g1=sub_samples.v.λ_g1[i], λ_g2=sub_samples.v.λ_g2[i],
                                     K_g=sub_samples.v.K_g[i], λ_q=sub_samples.v.λ_q[i], 
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
    
    plot(p1, p2, size=plot_size, xlabel="Bin number")
end
