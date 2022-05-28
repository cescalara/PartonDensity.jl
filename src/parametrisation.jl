using SpecialFunctions
using Parameters, Distributions
using Plots, Random
const sf = SpecialFunctions

export DirichletPDFParams, ValencePDFParams, UValencePDFParams
export DIRICHLET_TYPE, VALENCE_TYPE, UVALENCE_TYPE
export get_scaled_θ
export plot_input_pdfs, int_xtotx, xtotx
export get_input_pdf_func
export input_pdf_map

VALENCE_TYPE = 1
DIRICHLET_TYPE = 2
UVALENCE_TYPE = 3

"""
    Abstract type for any PDF parametrisation.
"""
abstract type AbstractPDFParams end

"""
    struct DirichletPDFParams <: AbstractPDFParams

Full Dirichlet specification of input PDF parameters.
"""
@with_kw struct DirichletPDFParams <: AbstractPDFParams
    param_type::Integer = DIRICHLET_TYPE
    seed::Integer = 0
    weights::Vector{Float64} = [1, 1, 1, 1, 1, 1, 1, 1, 1]
    θ::Vector{Float64} = rand(MersenneTwister(seed), Dirichlet(weights))
    K_u::Float64
    λ_u::Float64 = (θ[1] * (K_u + 1)) / (2 - θ[1])
    K_d::Float64
    λ_d::Float64 = (θ[2] * (K_d + 1)) / (1 - θ[2])
    λ_g1::Float64
    λ_g2::Float64
    K_g::Float64
    λ_q::Float64   
end

"""
    struct UValencePDFParams <: AbstractPDFParams

Full Dirichlet specification of input PDF parameters.
Includes 2 U-valence components.
"""
@with_kw struct UValencePDFParams <: AbstractPDFParams
    param_type::Integer = UVALENCE_TYPE
    seed::Integer = 0
    weights::Vector{Float64} = [1, 1, 1, 1, 1, 1, 1, 1, 1, 1]
    n_weights::Vector{Float64} = [1, 1]
    θ::Vector{Float64} = rand(MersenneTwister(seed), Dirichlet(weights))
    n::Vector{Float64} = 2 * rand(MersenneTwister(seed), Dirichlet(n_weights))
    K_u::Float64
    λ_u1::Float64 = (θ[1] * (K_u + 1)) / (n[1] - θ[1])
    λ_u2::Float64 = (θ[2] * (K_u + 1)) / (n[2] - θ[2])
    K_d::Float64
    λ_d::Float64 = (θ[3] * (K_d + 1)) / (1 - θ[3])
    λ_g1::Float64
    λ_g2::Float64
    K_g::Float64
    λ_q::Float64
end

"""
    struct ValencePDFParams <: AbstractPDFParams

Valence specification of input PDF parameters.
"""
@with_kw struct ValencePDFParams <: AbstractPDFParams
    param_type::Integer = VALENCE_TYPE
    λ_u::Float64
    K_u::Float64
    λ_d::Float64
    K_d::Float64
    λ_g1::Float64
    λ_g2::Float64
    K_g::Float64
    λ_q::Float64
    seed::Integer = 0
    weights::Vector{Float64} = [1, 1, 1, 1, 1, 1, 1]
    θ::Vector{Float64} = get_dirichlet_samples(λ_u, K_u, λ_d, K_d, seed, weights)
end


"""
    get_scaled_θ(λ_u, K_u, λ_d, K_d, θ_tmp)

Given a set of Dirichlet samples, θ_tmp, scale 
according to the valence params. Relevant for 
`ValencePDFParams`.
"""
function get_scaled_θ(λ_u::Float64, K_u::Float64, λ_d::Float64, K_d::Float64,
                      θ_tmp::Vector{Float64})

    
    A_u = 2 / sf.beta(λ_u, K_u + 1)
    A_d = 1 / sf.beta(λ_d, K_d + 1)

    I_u = A_u * sf.beta(λ_u + 1, K_u + 1)
    I_d = A_d * sf.beta(λ_d + 1, K_d + 1)

    remaining = 1 - (I_u + I_d)

    return θ_tmp * remaining
  
end


"""
    get_dirichlet_samples(λ_u, K_u, λ_d, K_d, seed, weights)

Given valance shape parameters and weights, get samples θ
that are correctly scaled. Relevant for `ValencePDFParams` 
"""
function get_dirichlet_samples(λ_u::Float64, K_u::Float64, λ_d::Float64,
                               K_d::Float64, seed::Integer, weights::Vector{Float64})

    Random.seed!(seed);
    θ_tmp = rand(Dirichlet(weights))

    θ = get_scaled_θ(λ_u, K_u, λ_d, K_d, θ_tmp)
    
    return θ
    
end


"""
    x_uv_x(x, λ_u, K_u)

Momentum density of u valence component.
Beta function 
    A_u x^λ_u (1-x)^K_u
A_u is set by λ_u and K_u.
"""
function x_uv_x(x::Float64, λ_u::Float64, K_u::Float64)
    
    A_u = 2 / sf.beta(λ_u, K_u + 1)
    
    return A_u * x^λ_u * (1 - x)^K_u
    
end

"""
    x_uv_x(x, λ_u1, λ_u2, K_u, n)

Momentum density of u valence component.
Beta function 
    A_u1 x^λ_u1 (1-x)^K_u
    A_u2 x^λ_u2 (1-x)^K_u
    
A_u is set by λ_u, K_u and n.

Alternative form for 2 components in u valence.
"""
function x_uv_x(x::Float64, λ_u1::Float64, λ_u2::Float64, K_u::Float64, n::Vector{Float64})

    A_u1 = n[1] / sf.beta(λ_u1, K_u + 1)
    A_u2 = n[2] / sf.beta(λ_u2, K_u + 1)

    comp1 = A_u1 * x^λ_u1 * (1 - x)^K_u
    comp2 = A_u2 * x^λ_u2 * (1 - x)^K_u
    
    return comp1 + comp2
    
end

"""
    x_dv_x(x, λ_d, K_d)

Momentum density of d valence component.
Beta function
    A_d x^λ_d (1 - x)^K_d
A_d is set by λ_d and K_d.
"""
function x_dv_x(x::Float64, λ_d::Float64, K_d::Float64)

    A_d = 1 / sf.beta(λ_d, K_d + 1)

    return A_d * x^λ_d * (1 - x)^K_d
    
end

"""
        x_g_x(x, λ_g1, λ_g2, K_g, w1, w2)

Momentum density of gluon component.
    A_g1 x^λ_g1 (1 - x)^K_g + A_g2 x^λ_g2
Amplitudes are set by weights `w1` and `w2`.
"""
function x_g_x(x::Float64, λ_g1::Float64, λ_g2::Float64, K_g::Float64,
               w1::Float64, w2::Float64)
    
    A_g1 = w1 / sf.beta(λ_g1 + 1, K_g + 1)
    
    x_g1_x = A_g1 * x^λ_g1 * (1 - x)^K_g

    A_g2 = w2 / sf.beta(λ_g2 + 1, 5 + 1)
    
    x_g2_x = A_g2 * x^λ_g2 * (1 - x)^5

    return x_g1_x + x_g2_x
    
end

"""
    x_q_x(x, λ_q, w)

Momentum density of non-valence
quark component.
    A_q x^λ_q
Amplitude is set by weight `w`
"""
function x_q_x(x::Float64, λ_q::Float64, w::Float64)
    
    A_q = (w / 2) / sf.beta(λ_q + 1, 5 + 1)
    
    return A_q * x^λ_q * (1 - x)^5
    
end


"""
    x_total_x(x, λ_u, K_u, λ_d, K_d, λ_g1, λ_g2, K_g, λ_q, θ)

Total momentum density.
"""
function xtotx(x::Float64, λ_u::Float64, K_u::Float64, λ_d::Float64, K_d::Float64,
               λ_g1::Float64, λ_g2::Float64, K_g::Float64, λ_q::Float64, θ::Vector{Float64})

    xuvx = x_uv_x(x, λ_u, K_u)

    xdvx = x_dv_x(x, λ_d, K_d)

    # For ValencePDFParams
    if length(θ) == 7
        
        xgx = x_g_x(x, λ_g1, λ_g2, K_g, θ[1], θ[2])

        xubarx = x_q_x(x, λ_q, θ[3])
        
        xdbarx = x_q_x(x, λ_q, θ[4])
        
        xsx = x_q_x(x, λ_q, θ[5])
        
        xcx = x_q_x(x, λ_q, θ[6])
    
        xbx = x_q_x(x, λ_q, θ[7])

    # For DirichletPDFParams
    elseif length(θ) == 9

        xgx = x_g_x(x, λ_g1, λ_g2, K_g, θ[3], θ[4])

        xubarx = x_q_x(x, λ_q, θ[5])
        
        xdbarx = x_q_x(x, λ_q, θ[6])
        
        xsx = x_q_x(x, λ_q, θ[7])
        
        xcx = x_q_x(x, λ_q, θ[8])
    
        xbx = x_q_x(x, λ_q, θ[9])

    else

        @error "length(θ) must be 7, 9 or 10 for a valid PDF parametrisation"

    end 
    
    return xuvx + xdvx + xgx + xubarx + xdbarx + xsx + xcx + xbx
    
end

"""
    x_total_x(x, λ_u1, λ_u2, K_u, λ_d, K_d, λ_g1, λ_g2, K_g, λ_q, θ)

Total momentum density.
"""
function xtotx(x::Float64, λ_u1::Float64, λ_u2::Float64, K_u::Float64, λ_d::Float64, K_d::Float64,
               λ_g1::Float64, λ_g2::Float64, K_g::Float64, λ_q::Float64, θ::Vector{Float64}, n::Vector{Float64})

    xuvx = x_uv_x(x, λ_u1, λ_u2, K_u, n)

    xdvx = x_dv_x(x, λ_d, K_d)

    xgx = x_g_x(x, λ_g1, λ_g2, K_g, θ[4], θ[5])
    
    xubarx = x_q_x(x, λ_q, θ[6])
    
    xdbarx = x_q_x(x, λ_q, θ[7])
        
    xsx = x_q_x(x, λ_q, θ[8])
    
    xcx = x_q_x(x, λ_q, θ[9])
    
    xbx = x_q_x(x, λ_q, θ[10])
    
    return xuvx + xdvx + xgx + xubarx + xdbarx + xsx + xcx + xbx
    
end


"""
    x_total_x(x, pdf_params)

Total momentum density.
"""
function xtotx(x::Float64, pdf_params::AbstractPDFParams)

    pdf = pdf_params

    if pdf.param_type == UVALENCE_TYPE

        return xtotx(x, pdf.λ_u1, pdf.λ_u2, pdf.K_u, pdf.λ_d, pdf.K_d, pdf.λ_g1,
                     pdf.λ_g2, pdf.K_g, pdf.λ_q, pdf.θ, pdf.n)

        
    else

        return xtotx(x, pdf.λ_u, pdf.K_u, pdf.λ_d, pdf.K_d, pdf.λ_g1,
                     pdf.λ_g2, pdf.K_g, pdf.λ_q, pdf.θ)

    end
        
end

"""
    int_xtotx(λ_u, K_u, λ_d, K_d, λ_g1, λ_g2, K_g, λ_q, θ)

Total integrated momentum density. Should equal 1.
"""
function int_xtotx(λ_u::Float64, K_u::Float64, λ_d::Float64, K_d::Float64,
                   λ_g1::Float64, λ_g2::Float64, K_g::Float64, λ_q::Float64, θ::Array{Float64})

    A_u = 2 / sf.beta(λ_u, K_u + 1)
    A_d = 1 / sf.beta(λ_d, K_d + 1)

    # ValencePDFParams
    if length(θ) == 7
        
        A_g1 = θ[1] / sf.beta(λ_g1 + 1, K_g + 1)
        A_g2 = θ[2] / sf.beta(λ_g2 + 1, 5 + 1)
    
        A_ubar = (θ[3] / 2) / sf.beta(λ_q + 1, 5 + 1)
        A_dbar = (θ[4] / 2) / sf.beta(λ_q + 1, 5 + 1)

        A_s = (θ[5] / 2) / sf.beta(λ_q + 1, 5 + 1)
        A_c = (θ[6] / 2) / sf.beta(λ_q + 1, 5 + 1)
        A_b = (θ[7] / 2) / sf.beta(λ_q + 1, 5 + 1)

    # DirichletPDFParams
    elseif length(θ) == 9

        A_g1 = θ[3] / sf.beta(λ_g1 + 1, K_g + 1)
        A_g2 = θ[4] / sf.beta(λ_g2 + 1, 5 + 1)
    
    
        A_ubar = (θ[5] / 2) / sf.beta(λ_q + 1, 5 + 1)
        A_dbar = (θ[6] / 2) / sf.beta(λ_q + 1, 5 + 1)

        A_s = (θ[7] / 2) / sf.beta(λ_q + 1, 5 + 1)
        A_c = (θ[8] / 2) / sf.beta(λ_q + 1, 5 + 1)
        A_b = (θ[9] / 2) / sf.beta(λ_q + 1, 5 + 1)
    
    else
        
         @error "length(θ) must be 7 or 9 for a valid PDF parametrisation"

    end       

    I_u = A_u * sf.beta(λ_u + 1, K_u + 1)
    I_d = A_d * sf.beta(λ_d + 1, K_d + 1)
    I_g = A_g1 * sf.beta(λ_g1 + 1, K_g + 1) + A_g2 * sf.beta(λ_g2 + 1, 5 + 1)
    I_q = 2 * (A_ubar + A_dbar + A_s + A_c + A_b) * sf.beta(λ_q + 1, 5 + 1)
    
    I_tot = I_u + I_d + I_g + I_q
    
    I_tot
end

"""
    int_xtotx(λ_u1, λ_u2, K_u, λ_d, K_d, λ_g1, λ_g2, K_g, λ_q, θ, n)

Total integrated momentum density. Should equal 1.
"""
function int_xtotx(λ_u1::Float64, λ_u2::Float64, K_u::Float64, λ_d::Float64, K_d::Float64,
                   λ_g1::Float64, λ_g2::Float64, K_g::Float64, λ_q::Float64, θ::Array{Float64},
                   n::Vector{Float64})

    A_u1 = n[1] / sf.beta(λ_u1, K_u + 1)
    A_u2 = n[2] / sf.beta(λ_u2, K_u + 1)
    A_d = 1 / sf.beta(λ_d, K_d + 1)

    A_g1 = θ[4] / sf.beta(λ_g1 + 1, K_g + 1)
    A_g2 = θ[5] / sf.beta(λ_g2 + 1, 5 + 1)
    
    
    A_ubar = (θ[6] / 2) / sf.beta(λ_q + 1, 5 + 1)
    A_dbar = (θ[7] / 2) / sf.beta(λ_q + 1, 5 + 1)

    A_s = (θ[8] / 2) / sf.beta(λ_q + 1, 5 + 1)
    A_c = (θ[9] / 2) / sf.beta(λ_q + 1, 5 + 1)
    A_b = (θ[10] / 2) / sf.beta(λ_q + 1, 5 + 1)
    
    I_u = A_u1 * sf.beta(λ_u1 + 1, K_u + 1) + A_u2 * sf.beta(λ_u2 + 1, K_u + 1)
    I_d = A_d * sf.beta(λ_d + 1, K_d + 1)
    I_g = A_g1 * sf.beta(λ_g1 + 1, K_g + 1) + A_g2 * sf.beta(λ_g2 + 1, 5 + 1)
    I_q = 2 * (A_ubar + A_dbar + A_s + A_c + A_b) * sf.beta(λ_q + 1, 5 + 1)
    
    I_tot = I_u + I_d + I_g + I_q
    
    I_tot
end


"""
    int_xtotx(pdf_params)

Total integrated momentum density. Should equal 1.
"""
function int_xtotx(pdf_params::AbstractPDFParams)

    pdf = pdf_params

    if pdf.param_type == UVALENCE_TYPE

        result = int_xtotx(pdf.λ_u1, pdf.λ_u2, pdf.K_u, pdf.λ_d, pdf.K_d,
                           pdf.λ_g1, pdf.λ_g2, pdf.K_g, pdf.λ_q, pdf.θ, pdf.n)

    else

        result = int_xtotx(pdf.λ_u, pdf.K_u, pdf.λ_d, pdf.K_d,
                           pdf.λ_g1, pdf.λ_g2, pdf.K_g, pdf.λ_q, pdf.θ)

    end
        
    return result
end

"""
    plot_input_pdfs(pdf_params, xmin, xmax, nx)

Plot the input PDFs defined by hyper_params over 
the given x range.   
"""
function plot_input_pdfs end


function plot_input_pdfs(pdf_params::ValencePDFParams, xmin::Float64=1.0e-2,
                         xmax::Float64=1.0, nx::Integer=1000)

    x_grid = range(xmin, stop=xmax, length=nx)
    pdf = pdf_params

    p = plot(x_grid, [x_uv_x(x, pdf.λ_u, pdf.K_u) for x in x_grid], label="x uv(x)", lw=3)
    p = plot!(x_grid, [x_dv_x(x, pdf.λ_d, pdf.K_d) for x in x_grid], label="x dv(x)", lw=3)

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


function plot_input_pdfs(pdf_params::DirichletPDFParams, xmin::Float64=1.0e-2,
                         xmax::Float64=1.0, nx::Integer=1000)

    x_grid = range(xmin, stop=xmax, length=nx)
    pdf = pdf_params

    p = plot(x_grid, [x_uv_x(x, pdf.λ_u, pdf.K_u) for x in x_grid], label="x uv(x)", lw=3)
    p = plot!(x_grid, [x_dv_x(x, pdf.λ_d, pdf.K_d) for x in x_grid], label="x dv(x)", lw=3)

    p = plot!(x_grid, [x_g_x(x, pdf.λ_g1, pdf.λ_g2, pdf.K_g, pdf.θ[3], pdf.θ[4]) 
                       for x in x_grid], label="x g(x)", lw=3)
    p = plot!(x_grid, [x_q_x(x, pdf.λ_q, pdf.θ[5]) for x in x_grid], label="x ubar(x)", lw=3)
    p = plot!(x_grid, [x_q_x(x, pdf.λ_q, pdf.θ[6]) for x in x_grid], label="x dbar(x)", lw=3)
    p = plot!(x_grid, [x_q_x(x, pdf.λ_q, pdf.θ[7]) for x in x_grid], label="x s(x)", lw=3)
    p = plot!(x_grid, [x_q_x(x, pdf.λ_q, pdf.θ[8]) for x in x_grid], label="x c(x)", lw=3)
    p = plot!(x_grid, [x_q_x(x, pdf.λ_q, pdf.θ[9]) for x in x_grid], label="x b(x)", lw=3)

    p = plot!(xlabel="x")
    p = plot!(xaxis=:log, yaxis=:log, legend=:outertopright)
    p = ylims!(1e-8, 30)

    return p
end


function plot_input_pdfs(pdf_params::UValencePDFParams, xmin::Float64=1.0e-2,
                         xmax::Float64=1.0, nx::Integer=1000)

    x_grid = range(xmin, stop=xmax, length=nx)
    pdf = pdf_params

    p = plot(x_grid, [x_uv_x(x, pdf.λ_u1, pdf.λ_u2, pdf.K_u, pdf.n) for x in x_grid], label="x uv(x)", lw=3)
    p = plot!(x_grid, [x_dv_x(x, pdf.λ_d, pdf.K_d) for x in x_grid], label="x dv(x)", lw=3)

    p = plot!(x_grid, [x_g_x(x, pdf.λ_g1, pdf.λ_g2, pdf.K_g, pdf.θ[4], pdf.θ[5]) 
                       for x in x_grid], label="x g(x)", lw=3)
    p = plot!(x_grid, [x_q_x(x, pdf.λ_q, pdf.θ[6]) for x in x_grid], label="x ubar(x)", lw=3)
    p = plot!(x_grid, [x_q_x(x, pdf.λ_q, pdf.θ[7]) for x in x_grid], label="x dbar(x)", lw=3)
    p = plot!(x_grid, [x_q_x(x, pdf.λ_q, pdf.θ[8]) for x in x_grid], label="x s(x)", lw=3)
    p = plot!(x_grid, [x_q_x(x, pdf.λ_q, pdf.θ[9]) for x in x_grid], label="x c(x)", lw=3)
    p = plot!(x_grid, [x_q_x(x, pdf.λ_q, pdf.θ[10]) for x in x_grid], label="x b(x)", lw=3)

    p = plot!(xlabel="x")
    p = plot!(xaxis=:log, yaxis=:log, legend=:outertopright)
    p = ylims!(1e-8, 30)

    return p
end

    
"""
    get_input_pdf_func(pdf_params)

Get the function to input the PDFs into QCDNUM.
"""
function get_input_pdf_func end


function get_input_pdf_func(pdf_params::ValencePDFParams)::Function

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
            f = x_uv_x(x, pdf.λ_u, pdf.K_u)
        end

        # d valence
        if (i == 2)
            f = x_dv_x(x, pdf.λ_d, pdf.K_d)
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


function get_input_pdf_func(pdf_params::DirichletPDFParams)::Function

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
            f = x_uv_x(x, pdf.λ_u, pdf.K_u)
        end

        # d valence
        if (i == 2)
            f = x_dv_x(x, pdf.λ_d, pdf.K_d)
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


function get_input_pdf_func(pdf_params::UValencePDFParams)::Function

    pdf = pdf_params
    
    func = function _input_pdfs(i, x)::Float64
        i = i[]
        x = x[]

        f = 0.0

        # gluon
        if (i == 0)
            f = x_g_x(x, pdf.λ_g1, pdf.λ_g2, pdf.K_g, pdf.θ[4], pdf.θ[5]) 
        end

        # u valence
        if (i == 1)
            f = x_uv_x(x, pdf.λ_u1, pdf.λ_u2, pdf.K_u, pdf.n)
        end

        # d valence
        if (i == 2)
            f = x_dv_x(x, pdf.λ_d, pdf.K_d)
        end

        # ubar
        if (i == 3)
            f = x_q_x(x, pdf.λ_q, pdf.θ[6])
        end

        # dbar
        if (i == 4)
            f = x_q_x(x, pdf.λ_q, pdf.θ[7])
        end

        # s and sbar
        if (i == 5) || (i == 6)
            f = x_q_x(x, pdf.λ_q, pdf.θ[8])
        end

        # c and cbar
        if (i == 7) || (i == 8)
            f = x_q_x(x, pdf.λ_q, pdf.θ[9])
        end

        # d and dbar
        if (i == 9) || (i == 10)
            f = x_q_x(x, pdf.λ_q, pdf.θ[9])
        end

        return f
    end

    return func
end

    
"""
    input_pdf_map

The relevant mapping for use with QCDNUM and `get_input_pdf_func`.
"""
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
