using SpecialFunctions
using Parameters, Distributions
using Plots, Random
const sf = SpecialFunctions

export BernsteinPDFParams
export BERNSTEIN_TYPE
export get_scaled_θ, get_scaled_UD
export plot_input_pdfs, int_xtotx, xtotx
export get_input_pdf_func
export input_pdf_map


BERNSTEIN_TYPE = 3


abstract type AbstractPDFParams end


@with_kw struct BernsteinPDFParams <: AbstractPDFParams
    param_type::Integer = BERNSTEIN_TYPE
    bspoly_params::Vector{Vector{Int64}} = [[0,3],[0,4],[1,4],[0,5]]
    U_weights::Vector{Float64} = zeros(length(bspoly_params))
    D_weights::Vector{Float64} = zeros(length(bspoly_params))
    λ_g1::Float64
    λ_g2::Float64
    K_g::Float64
    λ_q::Float64
    seed::Integer = 0
    weights::Vector{Float64} = ones(7)
    U_list::Vector{Float64} = get_dirichlet_UD(U_weights, 2, seed, bspoly_params)
    D_list::Vector{Float64} = get_dirichlet_UD(D_weights, 1, seed, bspoly_params)
    θ::Vector{Float64} = get_dirichlet_samples(U_list, D_list, seed, weights, bspoly_params)
end


function get_scaled_UD(UD_tmp::Vector{Float64}, intres::Integer, bspoly_params::Vector{Vector{Int64}} = [[0,3],[0,4],[1,4],[0,5]])
    
    intres_tmp = 0
    
    for (UD, params) in zip(UD_tmp, bspoly_params)
        intres_tmp += UD / (params[2] + 1)
    end
    
    return UD_tmp * (intres / intres_tmp)

end


function get_dirichlet_UD(UD_weights::Vector{Float64}, intres::Integer, seed::Integer, 
                          bspoly_params::Vector{Vector{Int64}})
    
    Random.seed!(seed);
    UD_tmp = rand(Dirichlet(UD_weights))
    
    return get_scaled_UD(UD_tmp, intres, bspoly_params)
    
end


function xfx_int(UD_list::Vector{Float64}, bspoly_params::Vector{Vector{Int64}})
        
    I = 0
        
    for (UD, params) in zip(UD_list, bspoly_params)
        I += UD * (params[1] + 1) / ((params[2] + 1)*(params[2] + 2))
    end
        
    return I
        
end


function get_scaled_θ(U_list::Vector{Float64}, D_list::Vector{Float64}, θ_tmp::Vector{Float64},
                      bspoly_params::Vector{Vector{Int64}} = [[0,3],[0,4],[1,4],[0,5]])
    
    I_u = xfx_int(U_list, bspoly_params)
    I_d = xfx_int(D_list, bspoly_params)

    remaining = 1 - (I_u + I_d)

    return θ_tmp * remaining
  
end


function get_dirichlet_samples(U_list::Vector{Float64}, D_list::Vector{Float64}, seed::Integer,
                                    weights::Vector{Float64}, bspoly_params::Vector{Vector{Int64}})

    Random.seed!(seed);
    θ_tmp = rand(Dirichlet(weights))

    θ = get_scaled_θ(U_list, D_list, θ_tmp, bspoly_params)
    
    if length(θ) != 7
        
        @error("length(θ) must be 7!")
        
    end
    
    return θ
    
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


function xtotx(x::Float64, pdf_params::BernsteinPDFParams)

    pdf = pdf_params 
    
    xux = x_uv_x(x, pdf.U_list, pdf.bspoly_params)
    
    xdx = x_dv_x(x, pdf.D_list, pdf.bspoly_params)
    
    xgx = x_g_x(x, pdf.λ_g1, pdf.λ_g2, pdf.K_g, pdf.θ[1], pdf.θ[2])
    
    xubarx = x_q_x(x, pdf.λ_q, pdf.θ[3])
    
    xdbarx = x_q_x(x, pdf.λ_q, pdf.θ[4])
    
    xsx = x_q_x(x, pdf.λ_q, pdf.θ[5])
    
    xcx = x_q_x(x, pdf.λ_q, pdf.θ[6])
    
    xbx = x_q_x(x, pdf.λ_q, pdf.θ[7])

    return xux + xdx + xgx + xubarx + xdbarx + xsx + xcx + xbx
end


function int_xtotx(pdf_params::BernsteinPDFParams)

    pdf = pdf_params
    
    I_u = xfx_int(pdf.U_list, pdf.bspoly_params)
    I_d = xfx_int(pdf.D_list, pdf.bspoly_params)
    
    A_g1 = pdf.θ[1] / sf.beta(pdf.λ_g1 + 1, pdf.K_g + 1)
    A_g2 = pdf.θ[2] / sf.beta(pdf.λ_g2 + 1, 5 + 1)
    
    A_ubar = (pdf.θ[3] / 2) / sf.beta(pdf.λ_q + 1, 5 + 1)
    A_dbar = (pdf.θ[4] / 2) / sf.beta(pdf.λ_q + 1, 5 + 1)

    A_s = (pdf.θ[5] / 2) / sf.beta(pdf.λ_q + 1, 5 + 1)
    A_c = (pdf.θ[6] / 2) / sf.beta(pdf.λ_q + 1, 5 + 1)
    A_b = (pdf.θ[7] / 2) / sf.beta(pdf.λ_q + 1, 5 + 1)
    
    I_g = A_g1 * sf.beta(pdf.λ_g1 + 1, pdf.K_g + 1) + A_g2 * sf.beta(pdf.λ_g2 + 1, 5 + 1)
    I_q = 2 * (A_ubar + A_dbar + A_s + A_c + A_b) * sf.beta(pdf.λ_q + 1, 5 + 1)
    
    return I_u + I_d + I_g + I_q
end


function plot_input_pdfs(pdf_params::BernsteinPDFParams, xmin::Float64=1.0e-2,
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


function get_input_pdf_func(pdf_params::BernsteinPDFParams)::Function

    pdf = pdf_params
    
    func = function _input_pdfs(i, x)::Float64
        i = i[]
        x = x[]

        f = 0.0

        # gluon
        if (i == 0)
            f = x_g_x(x, pdf.λ_g1, pdf.λ_g2, pdf.K_g, pdf.θ[1], pdf.θ[2]) 
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
            f = x_q_x(x, pdf.λ_q, pdf.θ[3])
        end

        # dbar
        if (i == 4)
            f = x_q_x(x, pdf.λ_q, pdf.θ[4])
        end

        # s and sbar
        if (i == 5) || (i == 6)
            f = x_q_x(x, pdf.λ_q, pdf.θ[5])
        end

        # c and cbar
        if (i == 7) || (i == 8)
            f = x_q_x(x, pdf.λ_q, pdf.θ[6])
        end

        # d and dbar
        if (i == 9) || (i == 10)
            f = x_q_x(x, pdf.λ_q, pdf.θ[7])
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
    
