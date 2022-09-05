export DirichletPDFParams, DIRICHLET_TYPE

DIRICHLET_TYPE = 2

"""
    struct DirichletPDFParams <: AbstractPDFParams

Full Dirichlet specification of input PDF parameters.

Constructors:

* ```$(FUNCTIONNAME)(; fields...)```

Fields:
$(TYPEDFIELDS)
"""
@with_kw struct DirichletPDFParams <: AbstractPDFParams
    param_type::Integer = DIRICHLET_TYPE
    seed::Integer = 0
    weights::Vector{Float64} = ones(9)
    θ::Vector{Float64} = rand(MersenneTwister(seed), Dirichlet(weights))
    K_u::Float64
    λ_u::Float64 = (θ[1] * (K_u + 1)) / (2 - θ[1])
    K_d::Float64
    λ_d::Float64 = (θ[2] * (K_d + 1)) / (1 - θ[2])
    λ_g1::Float64
    λ_g2::Float64
    K_g::Float64
    λ_q::Float64
    K_q::Float64
end

function xtotx(x::Float64, pdf_params::DirichletPDFParams)

    pdf = pdf_params

    return xtotx_dirichlet(x, pdf.λ_u, pdf.K_u, pdf.λ_d, pdf.K_d, pdf.λ_g1,
        pdf.λ_g2, pdf.K_g, pdf.λ_q, pdf.K_q, pdf.θ)

end

function xtotx_dirichlet(x::Float64, λ_u::Float64, K_u::Float64, λ_d::Float64,
    K_d::Float64, λ_g1::Float64, λ_g2::Float64, K_g::Float64, λ_q::Float64,
    K_q::Float64, θ::Vector{Float64})

    xuvx = x_uv_x(x, λ_u, K_u)

    xdvx = x_dv_x(x, λ_d, K_d)

    xgx = x_g_x(x, λ_g1, λ_g2, K_g, K_q, θ[3], θ[4])

    xubarx = x_q_x(x, λ_q, K_q, θ[5])

    xdbarx = x_q_x(x, λ_q, K_q, θ[6])

    xsx = x_q_x(x, λ_q, K_q, θ[7])

    xcx = x_q_x(x, λ_q, K_q, θ[8])

    xbx = x_q_x(x, λ_q, K_q, θ[9])

    return xuvx + xdvx + xgx + xubarx + xdbarx + xsx + xcx + xbx

end

function int_xtotx(pdf_params::DirichletPDFParams)

    pdf = pdf_params

    result = int_xtotx_dirichlet(pdf.λ_u, pdf.K_u, pdf.λ_d, pdf.K_d,
        pdf.λ_g1, pdf.λ_g2, pdf.K_g, pdf.λ_q, pdf.K_q, pdf.θ)

    return result
end

function int_xtotx_dirichlet(λ_u::Float64, K_u::Float64, λ_d::Float64,
    K_d::Float64, λ_g1::Float64, λ_g2::Float64, K_g::Float64, λ_q::Float64,
    K_q::Float64, θ::Array{Float64})

    A_u = 2 / sf.beta(λ_u, K_u + 1)
    A_d = 1 / sf.beta(λ_d, K_d + 1)

    A_g1 = θ[3] / sf.beta(λ_g1 + 1, K_g + 1)
    A_g2 = θ[4] / sf.beta(λ_g2 + 1, K_q + 1)

    A_ubar = (θ[5] / 2) / sf.beta(λ_q + 1, K_q + 1)
    A_dbar = (θ[6] / 2) / sf.beta(λ_q + 1, K_q + 1)

    A_s = (θ[7] / 2) / sf.beta(λ_q + 1, K_q + 1)
    A_c = (θ[8] / 2) / sf.beta(λ_q + 1, K_q + 1)
    A_b = (θ[9] / 2) / sf.beta(λ_q + 1, K_q + 1)

    I_u = A_u * sf.beta(λ_u + 1, K_u + 1)
    I_d = A_d * sf.beta(λ_d + 1, K_d + 1)
    I_g = A_g1 * sf.beta(λ_g1 + 1, K_g + 1) + A_g2 * sf.beta(λ_g2 + 1, K_q + 1)
    I_q = 2 * (A_ubar + A_dbar + A_s + A_c + A_b) * sf.beta(λ_q + 1, K_q + 1)

    I_tot = I_u + I_d + I_g + I_q

    I_tot
end

function plot_input_pdfs(pdf_params::DirichletPDFParams; xmin::Float64=1.0e-2,
    xmax::Float64=1.0, nx::Integer=1000)

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

function get_input_pdf_func(pdf_params::DirichletPDFParams)::Function

    pdf = pdf_params

    func = function _input_pdfs(i, x)::Float64
        i = i[]
        x = x[]

        f = 0.0

        # gluon
        if (i == 0)
            f = x_g_x(x, pdf.λ_g1, pdf.λ_g2, pdf.K_g, pdf.K_q, pdf.θ[3], pdf.θ[4])
        end

        # u valence
        if (i == 1)
            f = x_uv_x(x, pdf.λ_u, pdf.K_u)
        end

        # d valence
        if (i == 2)
            f = x_dv_x(x, pdf.λ_d, pdf.K_d)
        end

        # ubar
        if (i == 3)
            f = x_q_x(x, pdf.λ_q, pdf.K_q, pdf.θ[5])
        end

        # dbar
        if (i == 4)
            f = x_q_x(x, pdf.λ_q, pdf.K_q, pdf.θ[6])
        end

        # s and sbar
        if (i == 5) || (i == 6)
            f = x_q_x(x, pdf.λ_q, pdf.K_q, pdf.θ[7])
        end

        # c and cbar
        if (i == 7) || (i == 8)
            f = x_q_x(x, pdf.λ_q, pdf.K_q, pdf.θ[8])
        end

        # d and dbar
        if (i == 9) || (i == 10)
            f = x_q_x(x, pdf.λ_q, pdf.K_q, pdf.θ[9])
        end

        return f
    end

    return func
end