using SpecialFunctions
using Parameters, Distributions
using Plots, Random
const sf = SpecialFunctions

export BernDirPDFParams
export BERNDIR_TYPE
export get_UD_list
export plot_input_pdfs, int_xtotx, xtotx
export get_input_pdf_func

BERNDIR_TYPE = 4

abstract type AbstractPDFParams end

@with_kw struct BernDirPDFParams <: AbstractPDFParams
    param_type::Integer = BERNDIR_TYPE
    seed::Integer = 0
    weights::Vector{Float64} = ones(9)
    θ::Vector{Float64} = rand(MersenneTwister(seed), Dirichlet(weights))
    bspoly_params::Vector{Vector{Int64}} = [[0,3], [0,4], [1,4]]
    initial_U::Vector{Float64}
    initial_D::Vector{Float64}
    U_list::Vector{Float64} = get_UD_list(θ[1], 2, initial_U, bspoly_params)
    D_list::Vector{Float64} = get_UD_list(θ[2], 1, initial_D, bspoly_params)
    λ_g1::Float64
    λ_g2::Float64
    K_g::Float64
    λ_q::Float64   
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
        UD_2_extended += (xbxdx_factors[x+2] - bxdx_factors[x+2] * (xbxdx_factors[1]/bxdx_factors[1])) * 
                         initial_UD[x]
    end
    
    UD_2 = (Δudv - intres * (xbxdx_factors[1]/bxdx_factors[1]) - UD_2_extended) / 
           (xbxdx_factors[2] - bxdx_factors[2] * (xbxdx_factors[1]/bxdx_factors[1]))
    
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


function x_g_x(x::Float64, λ_g1::Float64, λ_g2::Float64, K_g::Float64,
               w1::Float64, w2::Float64)
    
    A_g1 = w1 / sf.beta(λ_g1 + 1, K_g + 1)
    
    x_g1_x = A_g1 * x^λ_g1 * (1 - x)^K_g

    A_g2 = w2 / sf.beta(λ_g2 + 1, 5 + 1)
    
    x_g2_x = A_g2 * x^λ_g2 * (1 - x)^5

    return x_g1_x + x_g2_x
    
end


function x_q_x(x::Float64, λ_q::Float64, w::Float64)
    
    A_q = (w / 2) / sf.beta(λ_q + 1, 5 + 1)
    
    return A_q * x^λ_q * (1 - x)^5
    
end


function xtotx(x::Float64, pdf_params::BernDirPDFParams)

    pdf = pdf_params 
    
    xux = x_uv_x(x, pdf.U_list, pdf.bspoly_params)
    
    xdx = x_dv_x(x, pdf.D_list, pdf.bspoly_params)
    
    xgx = x_g_x(x, pdf.λ_g1, pdf.λ_g2, pdf.K_g, pdf.θ[3], pdf.θ[4])
    
    xubarx = x_q_x(x, pdf.λ_q, pdf.θ[5])
    
    xdbarx = x_q_x(x, pdf.λ_q, pdf.θ[6])
    
    xsx = x_q_x(x, pdf.λ_q, pdf.θ[7])
    
    xcx = x_q_x(x, pdf.λ_q, pdf.θ[8])
    
    xbx = x_q_x(x, pdf.λ_q, pdf.θ[9])

    return xux + xdx + xgx + xubarx + xdbarx + xsx + xcx + xbx
end


function int_xtotx(pdf_params::BernDirPDFParams)

    pdf = pdf_params
    
    I_u = xfx_int(pdf.U_list, pdf.bspoly_params)
    I_d = xfx_int(pdf.D_list, pdf.bspoly_params)
    
    A_g1 = pdf.θ[3] / sf.beta(pdf.λ_g1 + 1, pdf.K_g + 1)
    A_g2 = pdf.θ[4] / sf.beta(pdf.λ_g2 + 1, 5 + 1)
    
    A_ubar = (pdf.θ[5] / 2) / sf.beta(pdf.λ_q + 1, 5 + 1)
    A_dbar = (pdf.θ[6] / 2) / sf.beta(pdf.λ_q + 1, 5 + 1)

    A_s = (pdf.θ[7] / 2) / sf.beta(pdf.λ_q + 1, 5 + 1)
    A_c = (pdf.θ[8] / 2) / sf.beta(pdf.λ_q + 1, 5 + 1)
    A_b = (pdf.θ[9] / 2) / sf.beta(pdf.λ_q + 1, 5 + 1)
    
    I_g = A_g1 * sf.beta(pdf.λ_g1 + 1, pdf.K_g + 1) + A_g2 * sf.beta(pdf.λ_g2 + 1, 5 + 1)
    I_q = 2 * (A_ubar + A_dbar + A_s + A_c + A_b) * sf.beta(pdf.λ_q + 1, 5 + 1)
    
    return I_u + I_d + I_g + I_q
end


function plot_input_pdfs(pdf_params::BernDirPDFParams, xmin::Float64=1.0e-2,
                         xmax::Float64=1.0, nx::Integer=1000)

    x_grid = range(xmin, stop=xmax, length=nx)
    pdf = pdf_params

    p = plot(x_grid, [x_uv_x(x, pdf.U_list, pdf.bspoly_params) for x in x_grid], label="x uv(x)", lw=3)
    p = plot!(x_grid, [x_dv_x(x, pdf.D_list, pdf.bspoly_params) for x in x_grid], label="x dv(x)", lw=3)

    p = plot!(x_grid, [x_g_x(x, pdf.λ_g1, pdf.λ_g2, pdf.K_g, pdf.θ[1], pdf.θ[2]) 
                       for x in x_grid], label="x g(x)", lw=3)    
    p = plot!(x_grid, [x_q_x(x, pdf.λ_q, pdf.θ[3]) for x in x_grid], label="x ubar(x)", lw=3)
    p = plot!(x_grid, [x_q_x(x, pdf.λ_q, pdf.θ[4]) for x in x_grid], label="x dbar(x)", lw=3)
    p = plot!(x_grid, [x_q_x(x, pdf.λ_q, pdf.θ[5]) for x in x_grid], label="x s(x)", lw=3)
    p = plot!(x_grid, [x_q_x(x, pdf.λ_q, pdf.θ[6]) for x in x_grid], label="x c(x)", lw=3)
    p = plot!(x_grid, [x_q_x(x, pdf.λ_q, pdf.θ[7]) for x in x_grid], label="x b(x)", lw=3)
        
    p = plot!(xlabel="x")
    p = plot!(xaxis=:log, yaxis=:log, legend=:outertopright)
    p = ylims!(1e-8, 30)

    return p
    
end


function get_input_pdf_func end


function get_input_pdf_func(pdf_params::BernDirPDFParams)::Function

    pdf = pdf_params
    
    func = function _input_pdfs(i, x)::Float64
        i = i[]
        x = x[]

        f = 0.0

        # gluon
        if (i == 0)
            f = x_g_x(x, pdf.λ_g1, pdf.λ_g2, pdf.K_g, pdf.θ[3], pdf.θ[4]) 
        end

        # u valence
        if (i == 1)
            f = x_uv_x(x, pdf.U_list, pdf.bspoly_params)
        end

        # d valence
        if (i == 2)
            f = x_dv_x(x, pdf.D_list, pdf.bspoly_params)
        end

        # ubar
        if (i == 3)
            f = x_q_x(x, pdf.λ_q, pdf.θ[5])
        end

        # dbar
        if (i == 4)
            f = x_q_x(x, pdf.λ_q, pdf.θ[6])
        end

        # s and sbar
        if (i == 5) || (i == 6)
            f = x_q_x(x, pdf.λ_q, pdf.θ[7])
        end

        # c and cbar
        if (i == 7) || (i == 8)
            f = x_q_x(x, pdf.λ_q, pdf.θ[8])
        end

        # d and dbar
        if (i == 9) || (i == 10)
            f = x_q_x(x, pdf.λ_q, pdf.θ[9])
        end

        return f
    end

    return func
end


#                         tb  bb  cb  sb  ub  db  g   d   u   s   c   b   t
input_pdf_map = Float64.([0., 0., 0., 0.,-1., 0., 0., 0., 1., 0., 0., 0., 0., # 1 # U valence
                          0., 0., 0., 0., 0.,-1., 0., 1., 0., 0., 0., 0., 0., # 2 # D valence
                          0., 0., 0., 0., 1., 0., 0., 0., 0., 0., 0., 0., 0., # 3 # u sea
                          0., 0., 0., 0., 0., 1., 0., 0., 0., 0., 0., 0., 0., # 4 # d sea
                          0., 0., 0., 0., 0., 0., 0., 0., 0., 1., 0., 0., 0., # 5 # s
                          0., 0., 0., 1., 0., 0., 0., 0., 0., 0., 0., 0., 0., # 6 # sbar
                          0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 1., 0., 0., # 7 # c
                          0., 0., 1., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., # 8 # cbar
                          0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 1., 0., # 9 # b
                          0., 1., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., # 10 # bbar
                          0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., # 11
                          0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.]); # 12
    
