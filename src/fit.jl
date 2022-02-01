export plot_model_space, plot_data_space

"""
    plot_model_space(pdf_params, samples)

Compare truth and posterior samples in the model space.
"""
function plot_model_space(pdf_params::PDFParameters, samples, xmin::Float64=1e-3,
                          xmax::Float64=1.0, nx::Integer=50, nsamples::Integer=200,
                          truth_color=:black, sample_color=:skyblue3)
    x_grid = range(xmin, stop=xmax, length=nx)

    p = plot(x_grid, [xtotx(x, pdf_params) for x in x_grid], color=truth_color,
             lw=3, label="Truth")

    for i in 1:nsamples

        pdf_params_i = PDFParameters(λ_u=samples.v.λ_u[i], K_u=samples.v.K_u[i], 
                                     λ_d=samples.v.λ_d[i], K_d=samples.v.K_d[i],
                                     λ_g1=samples.v.λ_g1[i], λ_g2=samples.v.λ_g2[i],
                                     K_g=samples.v.K_g[i], λ_q=samples.v.λ_q[i], 
                                     θ=Vector(samples.v.θ[i]))
        p = plot!(x_grid, [xtotx(x, pdf_params_i) for x in x_grid], color=sample_color, lw=3,
                  alpha=0.01, label="")
        
    end

    p = plot!(xaxis=:log, yaxis=:log)
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
                         splint_params::SPLINTParameters, quark_coeffs::QuarkCoefficients,
                         ep_color=:firebrick, em_color=:teal, nsamples::Integer=100)

    forward_model_init(qcdnum_grid, qcdnum_params, splint_params)

    nbins = sim_data["nbins"]

    p1 = scatter(1:nbins, sim_data["counts_obs_ep"], label="Observed counts (eP)", 
                 color=ep_color, lw=3)
    p2 = scatter(1:nbins, sim_data["counts_obs_ep"], label="Observed counts (eM)", 
                 color=em_color, lw=3)

    for i in 1:nsamples
        
        counts_obs_ep_i = zeros(UInt64, nbins)
        counts_obs_em_i = zeros(UInt64, nbins)

        pdf_params_i = PDFParameters(λ_u=samples.v.λ_u[i], K_u=samples.v.K_u[i], 
                                     λ_d=samples.v.λ_d[i], K_d=samples.v.K_d[i],
                                     λ_g1=samples.v.λ_g1[i], λ_g2=samples.v.λ_g2[i],
                                     K_g=samples.v.K_g[i], λ_q=samples.v.λ_q[i], 
                                     θ=Vector(samples.v.θ[i]))
        
        counts_pred_ep_i, counts_pred_em_i = forward_model(pdf_params_i, qcdnum_params, 
                                                           splint_params, quark_coeffs)
        
        for i in 1:nbins
            
            counts_obs_ep_i[i] = rand(Poisson(counts_pred_ep_i[i]))
            counts_obs_em_i[i] = rand(Poisson(counts_pred_em_i[i]))
            
        end
        
        p1 = scatter!(p1, 1:nbins, counts_obs_ep_i, label="", color=ep_color, 
                      lw=3, alpha=0.01)
        p2 = scatter!(p2, 1:nbins, counts_obs_em_i, label="", color=em_color, 
                      lw=3, alpha=0.01)
    
    end
    
    p1 = plot!(p1, xlabel="Bin number")
    p2 = plot!(p2, xlabel="Bin number")

    plot(p1, p2)
end
