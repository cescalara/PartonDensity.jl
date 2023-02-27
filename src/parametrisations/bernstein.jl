export BernsteinPDFParams, BERNSTEIN_TYPE
export BernsteinDirichletPDFParams, BERNSTEIN_DIRICHLET_TYPE

BERNSTEIN_TYPE = 3
BERNSTEIN_DIRICHLET_TYPE = 4


"""
    struct BernsteinPDFParams <: AbstractPDFParams

Bernstein polynomial specification of input PDF parameters.
Valence-style treatment. 

Constructors:

* ```$(FUNCTIONNAME)(; fields...)```

Fields:
$(TYPEDFIELDS)
"""
@with_kw struct BernsteinPDFParams <: AbstractPDFParams
    param_type::Int = BERNSTEIN_TYPE
    bspoly_params::Vector{Vector{Int64}} = [[0, 3], [0, 4], [1, 4], [0, 5]]
    bspoly_params_d::Vector{Vector{Int64}} = bspoly_params
    U_weights::Vector{Float64} = zeros(length(bspoly_params))
    D_weights::Vector{Float64} = zeros(length(bspoly_params))
    λ_g1::Float64
    λ_g2::Float64
    K_g::Float64
    λ_q::Float64
    K_q::Float64
    seed::Int = 0
    weights::Vector{Float64} = ones(7)
    U_list::Vector{Float64} = get_dirichlet_UD(U_weights, 2, seed, bspoly_params)
    D_list::Vector{Float64} = get_dirichlet_UD(D_weights, 1, seed, bspoly_params)
    θ::Vector{Float64} = get_dirichlet_samples(U_list, D_list, seed, weights, bspoly_params)
end

"""
    struct BernsteinDirichletPDFParams <: AbstractPDFParams

Bernstein polynomial specification of input PDF parameters.
Dirichlet-style treament. 

Constructors:

* ```$(FUNCTIONNAME)(; fields...)```

Fields:
$(TYPEDFIELDS)
"""
@with_kw struct BernsteinDirichletPDFParams <: AbstractPDFParams
    param_type::Int = BERNSTEIN_DIRICHLET_TYPE
    seed::Int = 0
    weights::Vector{Float64} = ones(9)
    θ::Vector{Float64} = rand(MersenneTwister(seed), Dirichlet(weights))
    bspoly_params::Vector{Vector{Int64}} = [[0, 3], [0, 4], [1, 4]]
    bspoly_params_d::Vector{Vector{Int64}} = bspoly_params
    initial_U::Vector{Float64} = ones(length(bspoly_params)-2)
    initial_D::Vector{Float64} = ones(length(bspoly_params_d)-2)
    U_list::Vector{Float64} = get_UD_list(θ[1], 2, initial_U, bspoly_params)
    D_list::Vector{Float64} = get_UD_list(θ[2], 1, initial_D, bspoly_params_d)
    λ_g1::Float64
    λ_g2::Float64
    K_g::Float64
    λ_q::Float64
    K_q::Float64
end

function get_scaled_UD(UD_tmp::Vector{Float64}, intres::Integer, bspoly_params::Vector{Vector{Int64}}=[[0, 3], [0, 4], [1, 4], [0, 5]])

    intres_tmp = 0

    for (UD, params) in zip(UD_tmp, bspoly_params)
        intres_tmp += UD / (params[2] + 1)
    end

    return UD_tmp * (intres / intres_tmp)

end

function get_dirichlet_UD(UD_weights::Vector{Float64}, intres::Integer, seed::Integer,
    bspoly_params::Vector{Vector{Int64}})

    Random.seed!(seed)
    UD_tmp = rand(Dirichlet(UD_weights))

    return get_scaled_UD(UD_tmp, intres, bspoly_params)

end

function xfx_int(UD_list::Vector{Float64}, bspoly_params::Vector{Vector{Int64}}=[[0, 3], [0, 4], [1, 4], [0, 5]])

    I = 0

    for (UD, params) in zip(UD_list, bspoly_params)
        I += UD * (params[1] + 1) / ((params[2] + 1) * (params[2] + 2))
    end

    return I

end

function get_scaled_θ(U_list::Vector{Float64}, D_list::Vector{Float64}, θ_tmp::Vector{Float64},
    bspoly_params::Vector{Vector{Int64}}=[[0, 3], [0, 4], [1, 4], [0, 5]])

    I_u = xfx_int(U_list, bspoly_params)
    I_d = xfx_int(D_list, bspoly_params)

    remaining = 1 - (I_u + I_d)

    return θ_tmp * remaining

end

function get_dirichlet_samples(U_list::Vector{Float64}, D_list::Vector{Float64}, seed::Integer,
    weights::Vector{Float64}, bspoly_params::Vector{Vector{Int64}})

    Random.seed!(seed)
    θ_tmp = rand(Dirichlet(weights))

    θ = get_scaled_θ(U_list, D_list, θ_tmp, bspoly_params)

    if length(θ) != 7

        @error("length(θ) must be 7!")

    end

    return θ

end

function get_bxdx(bspoly_params::Vector{Vector{Int64}})

    return [1 / (x[2] + 1) for x in bspoly_params]

end


function get_xbxdx(bspoly_params::Vector{Vector{Int64}})

    return [(x[1] + 1) / ((x[2] + 1) * (x[2] + 2)) for x in bspoly_params]

end


function get_UD_list(Δudv::Float64, intres::Int64, initial_UD::Vector{Float64},
    bspoly_params::Vector{Vector{Int64}})

    bxdx_factors = get_bxdx(bspoly_params)
    xbxdx_factors = get_xbxdx(bspoly_params)

    if length(initial_UD) != length(bspoly_params) - 2
        @error("The list for initial_UD must be 2 entries shorter than bspoly_params!")
    end

    UD_1_extended = 0
    UD_2_extended = 0

    for x in 1:length(initial_UD)
        UD_1_extended += bxdx_factors[x+2] * initial_UD[x]
        UD_2_extended += (xbxdx_factors[x+2] - bxdx_factors[x+2] * (xbxdx_factors[1] / bxdx_factors[1])) *
                         initial_UD[x]
    end

    UD_2 = (Δudv - intres * (xbxdx_factors[1] / bxdx_factors[1]) - UD_2_extended) /
           (xbxdx_factors[2] - bxdx_factors[2] * (xbxdx_factors[1] / bxdx_factors[1]))

    UD_1 = (intres - bxdx_factors[2] * UD_2 - UD_1_extended) / bxdx_factors[1]

    return append!([UD_1, UD_2], initial_UD)

end

function bs_poly(x::Float64, vn::Vector{Int64})

    v = vn[1]
    n = vn[2]

    return binomial(n, v) * x^v * (1 - x)^(n - v)

end

function x_uv_x(x::Float64, U_list::Vector{Float64}, bspoly_params::Vector{Vector{Int64}})

    res = 0

    for (U, params) in zip(U_list, bspoly_params)
        res += U * bs_poly(x, params)
    end

    return x * res

end

function x_dv_x(x::Float64, D_list::Vector{Float64}, bspoly_params::Vector{Vector{Int64}})

    res = 0

    for (D, params) in zip(D_list, bspoly_params)
        res += D * bs_poly(x, params)
    end

    return x * res

end

function xtotx(x::Float64, pdf_params::Union{BernsteinPDFParams,BernsteinDirichletPDFParams})

    pdf = pdf_params

    xux = x_uv_x(x, pdf.U_list, pdf.bspoly_params)

    xdx = x_dv_x(x, pdf.D_list, pdf.bspoly_params_d)

    xgx = x_g_x(x, pdf.λ_g1, pdf.λ_g2, pdf.K_g, pdf.K_q, pdf.θ[1], pdf.θ[2])

    xubarx = x_q_x(x, pdf.λ_q, pdf.K_q, pdf.θ[3])

    xdbarx = x_q_x(x, pdf.λ_q, pdf.K_q, pdf.θ[4])

    xsx = x_q_x(x, pdf.λ_q, pdf.K_q, pdf.θ[5])

    xcx = x_q_x(x, pdf.λ_q, pdf.K_q, pdf.θ[6])

    xbx = x_q_x(x, pdf.λ_q, pdf.K_q, pdf.θ[7])

    return xux + xdx + xgx + xubarx + xdbarx + xsx + xcx + xbx
end


function int_xtotx(pdf_params::Union{BernsteinPDFParams,BernsteinDirichletPDFParams})

    pdf = pdf_params

    I_u = xfx_int(pdf.U_list, pdf.bspoly_params)
    I_d = xfx_int(pdf.D_list, pdf.bspoly_params_d)

    A_g1 = pdf.θ[1] / sf.beta(pdf.λ_g1 + 1, pdf.K_g + 1)
    A_g2 = pdf.θ[2] / sf.beta(pdf.λ_g2 + 1, pdf.K_q + 1)

    A_ubar = (pdf.θ[3] / 2) / sf.beta(pdf.λ_q + 1, pdf.K_q + 1)
    A_dbar = (pdf.θ[4] / 2) / sf.beta(pdf.λ_q + 1, pdf.K_q + 1)

    A_s = (pdf.θ[5] / 2) / sf.beta(pdf.λ_q + 1, pdf.K_q + 1)
    A_c = (pdf.θ[6] / 2) / sf.beta(pdf.λ_q + 1, pdf.K_q + 1)
    A_b = (pdf.θ[7] / 2) / sf.beta(pdf.λ_q + 1, pdf.K_q + 1)

    I_g = A_g1 * sf.beta(pdf.λ_g1 + 1, pdf.K_g + 1) + A_g2 * sf.beta(pdf.λ_g2 + 1, pdf.K_q + 1)
    I_q = 2 * (A_ubar + A_dbar + A_s + A_c + A_b) * sf.beta(pdf.λ_q + 1, pdf.K_q + 1)

    return I_u + I_d + I_g + I_q
end

function plot_input_pdfs(pdf_params::Union{BernsteinPDFParams,BernsteinDirichletPDFParams};
    xmin::Float64=1.0e-2, xmax::Float64=1.0, nx::Integer=1000)

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

function get_input_pdf_func(pdf_params::Union{BernsteinPDFParams,BernsteinDirichletPDFParams})::Function

    pdf = pdf_params

    func = function _input_pdfs(i, x)::Float64
        i = i[]
        x = x[]

        f = 0.0

        # gluon
        if (i == 0)
            f = x_g_x(x, pdf.λ_g1, pdf.λ_g2, pdf.K_g, pdf.K_q, pdf.θ[1], pdf.θ[2])
        end

        # u valence
        if (i == 1)
            f = x_uv_x(x, pdf.U_list, pdf.bspoly_params)
        end

        # d valence
        if (i == 2)
            f = x_dv_x(x, pdf.D_list, pdf.bspoly_params_d)
        end

        # ubar
        if (i == 3)
            f = x_q_x(x, pdf.λ_q, pdf.K_q, pdf.θ[3])
        end

        # dbar
        if (i == 4)
            f = x_q_x(x, pdf.λ_q, pdf.K_q, pdf.θ[4])
        end

        # s and sbar
        if (i == 5) || (i == 6)
            f = x_q_x(x, pdf.λ_q, pdf.K_q, pdf.θ[5])
        end

        # c and cbar
        if (i == 7) || (i == 8)
            f = x_q_x(x, pdf.λ_q, pdf.K_q, pdf.θ[6])
        end

        # d and dbar
        if (i == 9) || (i == 10)
            f = x_q_x(x, pdf.λ_q, pdf.K_q, pdf.θ[7])
        end

        return f
    end

    return func
end
