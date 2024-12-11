# This file is a part of PartonDensity.jl, licensed under the MIT License (MIT).

module PartonDensityPlotsExt

using PartonDensity

using Plots

using PartonDensity: AbstractPDFParams
using QCDNUM

using PartonDensity: x_uv_x, x_dv_x, x_g_x, x_q_x, x_q_x, x_q_x, x_q_x, x_q_x


function PartonDensity.plot_input_pdfs(pdf_params::Union{BernsteinPDFParams,BernsteinDirichletPDFParams};
    xmin::Real=1.0e-2, xmax::Real=1.0, nx::Integer=1000)

    x_grid = range(xmin, stop=xmax, length=nx)
    pdf = pdf_params

    p = plot(x_grid, [x_uv_x(x, pdf.U_list, pdf.bspoly_params) for x in x_grid], label="x uv(x)", lw=3)
    p = plot!(x_grid, [x_dv_x(x, pdf.D_list, pdf.bspoly_params_d) for x in x_grid], label="x dv(x)", lw=3)

    p = plot!(x_grid, [x_g_x(x, pdf.λ_g1, pdf.λ_g2, pdf.K_g, pdf.K_q, pdf.θ[1], pdf.θ[2])
                       for x in x_grid], label="x g(x)", lw=3)
    p = plot!(x_grid, [x_q_x(x, pdf.λ_q, pdf.K_q, pdf.θ[3]) for x in x_grid], label="x ubar(x)", lw=3)
    p = plot!(x_grid, [x_q_x(x, pdf.λ_q, pdf.K_q, pdf.θ[4]) for x in x_grid], label="x dbar(x)", lw=3)
    p = plot!(x_grid, [x_q_x(x, pdf.λ_q, pdf.K_q, pdf.θ[5]) for x in x_grid], label="x s(x)", lw=3)
    p = plot!(x_grid, [x_q_x(x, pdf.λ_q, pdf.K_q, pdf.θ[6]) for x in x_grid], label="x c(x)", lw=3)
    p = plot!(x_grid, [x_q_x(x, pdf.λ_q, pdf.K_q, pdf.θ[7]) for x in x_grid], label="x b(x)", lw=3)

    p = plot!(xlabel="x")
    p = plot!(xaxis=:log, yaxis=:log, legend=:outertopright)
    p = ylims!(1e-8, 30)

    return p

end


function PartonDensity.plot_input_pdfs(pdf_params::ValencePDFParams; xmin::Float64=1.0e-2,
    xmax::Float64=1.0, nx::Integer=1000)

    x_grid = range(xmin, stop=xmax, length=nx)
    pdf = pdf_params

    p = plot(x_grid, [x_uv_x(x, pdf.λ_u, pdf.K_u) for x in x_grid], label="x uv(x)", lw=3)
    p = plot!(x_grid, [x_dv_x(x, pdf.λ_d, pdf.K_d) for x in x_grid], label="x dv(x)", lw=3)

    p = plot!(x_grid, [x_g_x(x, pdf.λ_g1, pdf.λ_g2, pdf.K_g, pdf.K_q, pdf.θ[1], pdf.θ[2])
                       for x in x_grid], label="x g(x)", lw=3)
    p = plot!(x_grid, [x_q_x(x, pdf.λ_q, pdf.K_q, pdf.θ[3]) for x in x_grid], label="x ubar(x)", lw=3)
    p = plot!(x_grid, [x_q_x(x, pdf.λ_q, pdf.K_q, pdf.θ[4]) for x in x_grid], label="x dbar(x)", lw=3)
    p = plot!(x_grid, [x_q_x(x, pdf.λ_q, pdf.K_q, pdf.θ[5]) for x in x_grid], label="x s(x)", lw=3)
    p = plot!(x_grid, [x_q_x(x, pdf.λ_q, pdf.K_q, pdf.θ[6]) for x in x_grid], label="x c(x)", lw=3)
    p = plot!(x_grid, [x_q_x(x, pdf.λ_q, pdf.K_q, pdf.θ[7]) for x in x_grid], label="x b(x)", lw=3)

    p = plot!(xlabel="x")
    p = plot!(xaxis=:log, yaxis=:log, legend=:outertopright)
    p = ylims!(1e-8, 30)

    return p
end


function PartonDensity.plot_input_pdfs(pdf_params::DirichletPDFParams; xmin::Real=1.0e-2,
    xmax::Real=1.0, nx::Integer=1000)

    x_grid = range(xmin, stop=xmax, length=nx)
    pdf = pdf_params

    p = plot(x_grid, [x_uv_x(x, pdf.λ_u, pdf.K_u) for x in x_grid], label="x uv(x)", lw=3)
    p = plot!(x_grid, [x_dv_x(x, pdf.λ_d, pdf.K_d) for x in x_grid], label="x dv(x)", lw=3)

    p = plot!(x_grid, [x_g_x(x, pdf.λ_g1, pdf.λ_g2, pdf.K_g, pdf.K_q, pdf.θ[3], pdf.θ[4])
                       for x in x_grid], label="x g(x)", lw=3)
    p = plot!(x_grid, [x_q_x(x, pdf.λ_q, pdf.K_q, pdf.θ[5]) for x in x_grid], label="x ubar(x)", lw=3)
    p = plot!(x_grid, [x_q_x(x, pdf.λ_q, pdf.K_q, pdf.θ[6]) for x in x_grid], label="x dbar(x)", lw=3)
    p = plot!(x_grid, [x_q_x(x, pdf.λ_q, pdf.K_q, pdf.θ[7]) for x in x_grid], label="x s(x)", lw=3)
    p = plot!(x_grid, [x_q_x(x, pdf.λ_q, pdf.K_q, pdf.θ[8]) for x in x_grid], label="x c(x)", lw=3)
    p = plot!(x_grid, [x_q_x(x, pdf.λ_q, pdf.K_q, pdf.θ[9]) for x in x_grid], label="x b(x)", lw=3)

    p = plot!(xlabel="x")
    p = plot!(xaxis=:log, yaxis=:log, legend=:outertopright)
    p = ylims!(1e-8, 30)

    return p
end


function PartonDensity.plot_model_space(pdf_params::AbstractPDFParams, samples; xmin::Float64=1e-3,
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


function plot_model_space_impl(x_grid::StepRangeLen{Float64}, pdf_params::ValencePDFParams, samples, p; color=:skyblue3)

    for i in eachindex(samples)

        θ_i = get_scaled_θ(samples.v.λ_u[i], samples.v.K_u[i], samples.v.λ_d[i],
            samples.v.K_d[i], Vector(samples.v.θ_tmp[i]))

        pdf_params_i = ValencePDFParams(λ_u=samples.v.λ_u[i], K_u=samples.v.K_u[i],
            λ_d=samples.v.λ_d[i], K_d=samples.v.K_d[i],
            λ_g1=samples.v.λ_g1[i], λ_g2=samples.v.λ_g2[i],
            K_g=samples.v.K_g[i], λ_q=samples.v.λ_q[i], K_q=samples.v.K_q[i],
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
            K_g=samples.v.K_g[i], λ_q=samples.v.λ_q[i], K_q=samples.v.K_q[i],
            θ=Vector(samples.v.θ[i]))
        p = plot!(x_grid, [xtotx(x, pdf_params_i) for x in x_grid], color=color, lw=3,
            alpha=0.01, label="")

    end

    return p
end


function PartonDensity.plot_data_space(pdf_params::AbstractPDFParams, sim_data::Dict{String,Any}, samples,
    qcdnum_params::QCDNUM.EvolutionParams,
    splint_params::QCDNUM.SPLINTParams, quark_coeffs::QuarkCoefficients, md::MetaData;
    ep_color=:firebrick, em_color=:teal, nsamples::Integer=100, plot_size=(1000, 500),
    args...)

    forward_model_init(qcdnum_params, splint_params)

    nbins = sim_data["nbins"]

    p1 = scatter(1:nbins, sim_data["counts_obs_ep"], label="Observed counts (eP)",
        color=ep_color, lw=3)
    p2 = scatter(1:nbins, sim_data["counts_obs_em"], label="Observed counts (eM)",
        color=em_color, lw=3)

    sub_samples = BAT.bat_sample(samples, BAT.OrderedResampling(nsamples=nsamples)).result

    p1, p2 = plot_data_space_impl(pdf_params, sub_samples, qcdnum_params,
        splint_params, quark_coeffs, md, p1, p2, nbins)

    plot(p1, p2, size=plot_size, xlabel="Bin number", bottom_margin=10Plots.mm, markerstrokewidth=0; args...)
end


function plot_data_space_impl(pdf_params::ValencePDFParams, samples, qcdnum_params::QCDNUM.EvolutionParams,
    splint_params::QCDNUM.SPLINTParams, quark_coeffs::QuarkCoefficients, md::MetaData,
    p1, p2, nbins::Integer; ep_color=:firebrick, em_color=:teal)

    for i in eachindex(samples)

        counts_obs_ep_i = zeros(UInt64, nbins)
        counts_obs_em_i = zeros(UInt64, nbins)

        θ_i = get_scaled_θ(samples.v.λ_u[i], samples.v.K_u[i], samples.v.λ_d[i],
            samples.v.K_d[i], Vector(samples.v.θ_tmp[i]))


        pdf_params_i = ValencePDFParams(λ_u=samples.v.λ_u[i], K_u=samples.v.K_u[i],
            λ_d=samples.v.λ_d[i], K_d=samples.v.K_d[i],
            λ_g1=samples.v.λ_g1[i], λ_g2=samples.v.λ_g2[i],
            K_g=samples.v.K_g[i], λ_q=samples.v.λ_q[i], K_q=samples.v.K_q[i],
            θ=θ_i)

        counts_pred_ep_i, counts_pred_em_i = forward_model(pdf_params_i, qcdnum_params, splint_params, quark_coeffs, md)

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


function plot_data_space_impl(pdf_params::DirichletPDFParams, samples, qcdnum_params::QCDNUM.EvolutionParams,
    splint_params::QCDNUM.SPLINTParams, quark_coeffs::QuarkCoefficients, md::MetaData,
    p1, p2, nbins::Integer; ep_color=:firebrick, em_color=:teal)

    for i in eachindex(samples)

        counts_obs_ep_i = zeros(UInt64, nbins)
        counts_obs_em_i = zeros(UInt64, nbins)

        pdf_params_i = DirichletPDFParams(K_u=samples.v.K_u[i], K_d=samples.v.K_d[i],
            λ_g1=samples.v.λ_g1[i], λ_g2=samples.v.λ_g2[i],
            K_g=samples.v.K_g[i], λ_q=samples.v.λ_q[i], K_q=samples.v.K_q[i],
            θ=Vector(samples.v.θ[i]))

        counts_pred_ep_i, counts_pred_em_i = forward_model(pdf_params_i, qcdnum_params, splint_params, quark_coeffs, md)

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

        p1 = scatter!(p1, 1:nbins, counts_obs_ep_i, label="", color=ep_color, lw=3, alpha=0.01)
        p2 = scatter!(p2, 1:nbins, counts_obs_em_i, label="", color=em_color, lw=3, alpha=0.01)

    end

    return p1, p2
end
########################

function plot_model_space_impl(x_grid::StepRangeLen{Float64}, pdf_params::BernsteinPDFParams, samples, p; color=:skyblue3)

    for i in eachindex(samples)

        U_list = get_scaled_UD(Vector(samples.v.U_weights[i]), 2)
        D_list = get_scaled_UD(Vector(samples.v.U_weights[i]), 1)

        θ_i = get_scaled_θ(U_list, D_list, Vector(samples.v.θ_tmp[i]))

        vec_bspp = Vector(samples.bspoly_params[i])
        bspoly_params = [[vec_bspp[Int(2 * i - 1)], vec_bspp[Int(2 * i)]] for i in 1:length(vec_bspp)/2]

        bspoly_params_d = 0

        try
            vec_bsppd = Vector(samples.bspoly_params_d[i])
            bspoly_params_d = [[vec_bsppd[Int(2 * i - 1)], vec_bsppd[Int(2 * i)]] for i in 1:length(vec_bsppd)/2]
        catch err
            bspoly_params_d = bspoly_params
        end

        pdf_params_i = BernsteinPDFParams(U_list=U_list, D_list=D_list,
            λ_g1=samples.v.λ_g1[i], λ_g2=samples.v.λ_g2[i],
            K_g=samples.v.K_g[i], λ_q=samples.v.λ_q[i], K_q=samples.v.K_q[i],
            θ=θ_i,
            bspoly_params=bspoly_params,
            bspoly_params_d=bspoly_params_d)

        p = plot!(x_grid, [xtotx(x, pdf_params_i) for x in x_grid], color=color, lw=3,
            alpha=0.01, label="")

    end

    return p
end


function plot_model_space_impl(x_grid::StepRangeLen{Float64}, pdf_params::BernsteinDirichletPDFParams, samples, p; color=:skyblue3)

    for i in eachindex(samples)

        vec_bspp = Vector(samples.bspoly_params[i])
        bspoly_params = [[vec_bspp[Int(2 * i - 1)], vec_bspp[Int(2 * i)]] for i in 1:length(vec_bspp)/2]

        bspoly_params_d = 0

        try
            vec_bsppd = Vector(samples.bspoly_params_d[i])
            bspoly_params_d = [[vec_bsppd[Int(2 * i - 1)], vec_bsppd[Int(2 * i)]] for i in 1:length(vec_bsppd)/2]
        catch err
            bspoly_params_d = bspoly_params
        end

        pdf_params_i = BernsteinDirichletPDFParams(initial_U=[samples.v.initial_U[i]], initial_D=[samples.v.initial_D[i]],
            λ_g1=samples.v.λ_g1[i], λ_g2=samples.v.λ_g2[i],
            K_g=samples.v.K_g[i], λ_q=samples.v.λ_q[i],
            θ=Vector(samples.v.θ[i]),
            bspoly_params=bspoly_params,
            bspoly_params_d=bspoly_params_d)

        p = plot!(x_grid, [xtotx(x, pdf_params_i) for x in x_grid], color=color, lw=3,
            alpha=0.01, label="")

    end

    return p
end

function plot_data_space_impl(pdf_params::BernsteinPDFParams, samples, qcdnum_params::QCDNUM.EvolutionParams,
    splint_params::QCDNUM.SPLINTParams, quark_coeffs::QuarkCoefficients,
    p1, p2, nbins::Integer; ep_color=:firebrick, em_color=:teal, md::MetaData)

    for i in eachindex(samples)

        counts_obs_ep_i = zeros(UInt64, nbins)
        counts_obs_em_i = zeros(UInt64, nbins)

        U_list = get_scaled_UD(Vector(samples.v.U_weights[i]), 2)
        D_list = get_scaled_UD(Vector(samples.v.U_weights[i]), 1)

        θ_i = get_scaled_θ(U_list, D_list, Vector(samples.v.θ_tmp[i]))

        vec_bspp = Vector(samples.bspoly_params[i])
        bspoly_params = [[vec_bspp[Int(2 * i - 1)], vec_bspp[Int(2 * i)]] for i in 1:length(vec_bspp)/2]

        bspoly_params_d = 0

        try
            vec_bsppd = Vector(samples.bspoly_params_d[i])
            bspoly_params_d = [[vec_bsppd[Int(2 * i - 1)], vec_bsppd[Int(2 * i)]] for i in 1:length(vec_bsppd)/2]
        catch err
            bspoly_params_d = bspoly_params
        end

        pdf_params_i = BernsteinPDFParams(U_list=U_list, D_list=D_list,
            λ_g1=samples.v.λ_g1[i], λ_g2=samples.v.λ_g2[i],
            K_g=samples.v.K_g[i], λ_q=samples.v.λ_q[i], K_q=samples.v.K_q[i],
            θ=θ_i, bspoly_params=bspoly_params,
            bspoly_params_d=bspoly_params_d)

        counts_pred_ep_i, counts_pred_em_i = forward_model(pdf_params_i, qcdnum_params, splint_params, quark_coeffs, md)

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


function plot_data_space_impl(pdf_params::BernsteinDirichletPDFParams, samples, qcdnum_params::QCDNUM.EvolutionParams,
    splint_params::QCDNUM.SPLINTParams, quark_coeffs::QuarkCoefficients,
    p1, p2, nbins::Integer; ep_color=:firebrick, em_color=:teal, md::MetaData)

    for i in eachindex(samples)

        counts_obs_ep_i = zeros(UInt64, nbins)
        counts_obs_em_i = zeros(UInt64, nbins)

        vec_bspp = Vector(samples.bspoly_params[i])
        bspoly_params = [[vec_bspp[Int(2 * i - 1)], vec_bspp[Int(2 * i)]] for i in 1:length(vec_bspp)/2]

        bspoly_params_d = 0

        try
            vec_bsppd = Vector(samples.bspoly_params_d[i])
            bspoly_params_d = [[vec_bsppd[Int(2 * i - 1)], vec_bsppd[Int(2 * i)]] for i in 1:length(vec_bsppd)/2]
        catch err
            bspoly_params_d = bspoly_params
        end

        pdf_params_i = BernsteinDirichletPDFParams(initial_U=[samples.v.initial_U[i]], initial_D=[samples.v.initial_D[i]],
            λ_g1=samples.v.λ_g1[i], λ_g2=samples.v.λ_g2[i],
            K_g=samples.v.K_g[i], λ_q=samples.v.λ_q[i],
            θ=Vector(samples.v.θ[i]),
            bspoly_params=bspoly_params,
            bspoly_params_d=bspoly_params_d)

        counts_pred_ep_i, counts_pred_em_i = forward_model(pdf_params_i, qcdnum_params, splint_params, quark_coeffs, md)

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








end # module PartonDensityPlotsExt
