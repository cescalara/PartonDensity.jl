using SpecialFunctions, Roots
using Plots
const sf = SpecialFunctions

"""
    x_uv_x(x, λ_u, w)

Momentum density of u valence component.
Beta function 
    A_u x^λ_u (1-x)^K_u
A_u and K_u are set by the passed weight, 
`w`. 
"""
function x_uv_x(x::Float64, λ_u::Float64, w::Float64)
    
    f(K_u) = w - 2 * (sf.beta(λ_u + 1, K_u + 1) / sf.beta(λ_u, K_u + 1))
    
    K_u = find_zero(f, (-1.1, 100))
    
    A_u = 2 / sf.beta(λ_u, K_u + 1)
    
    return A_u * x^λ_u * (1 - x)^K_u
    
end

"""
    x_dv_x(x, λ_d, w)

Momentum density of d valence component.
Beta function
    A_d x^λ_d (1 - x)^K_d
A_d and K_d are set by the passed weight,
`w`
"""
function x_dv_x(x::Float64, λ_d::Float64, w::Float64)

    f(K_d) = w - (sf.beta(λ_d + 1, K_d + 1) / sf.beta(λ_d, K_d + 1))
    
    K_d = find_zero(f, 1)

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
    x_total_x(x, λ_u, λ_d, λ_g1, λ_g2, K_g, λ_q, θ)

Total momentum density.
"""
function xtotx(x::Float64, λ_u::Float64, λ_d::Float64, λ_g1::Float64, λ_g2::Float64,
               K_g::Float64, λ_q::Float64, θ::Array{Float64})
    
    xuvx = x_uv_x(x, λ_u, θ[1])

    xdvx = x_dv_x(x, λ_d, θ[2])

    xgx = x_g_x(x, λ_g1, λ_g2, K_g, θ[3], θ[4])

    xubarx = x_q_x(x, λ_q, θ[5])
    
    xdbarx = x_q_x(x, λ_q, θ[6])

    xsx = x_q_x(x, λ_q, θ[7])
    
    xcx = x_q_x(x, λ_q, θ[8])
    
    xbx = x_q_x(x, λ_q, θ[9])
    
    return xuvx + xdvx + xgx + xubarx + xdbarx + xsx + xcx + xbx
    
end

"""
    int_xtotx(λ_u, λ_d, λ_g1, λ_g2, K_g, λ_q, θ)

Total integrated momentum density. Should equal 1.
"""
function int_xtotx(λ_u::Float64, λ_d::Float64, λ_g1::Float64, λ_g2::Float64,
                            K_g::Float64, λ_q::Float64, θ::Array{Float64})

    # Amplitudes
    fu(K_u) = θ[1] - 2 * (sf.beta(λ_u + 1, K_u + 1) / sf.beta(λ_u, K_u + 1))
    K_u = find_zero(fu, 1)
    A_u = 2 / sf.beta(λ_u, K_u + 1)

    fd(K_d) = θ[2] - (sf.beta(λ_d + 1, K_d + 1) / sf.beta(λ_d, K_d + 1))
    K_d = find_zero(fd, 1)
    A_d = 1 / sf.beta(λ_d, K_d + 1)

    A_g1 = θ[3] / sf.beta(λ_g1 + 1, K_g + 1)
    A_g2 = θ[4] / sf.beta(λ_g2 + 1, 5 + 1)
    
    
    A_ubar = (θ[5] / 2) / sf.beta(λ_q + 1, 5 + 1)

    A_dbar = (θ[6] / 2) / sf.beta(λ_q + 1, 5 + 1)

    A_s = (θ[7] / 2) / sf.beta(λ_q + 1, 5 + 1)

    A_c = (θ[8] / 2) / sf.beta(λ_q + 1, 5 + 1)

    A_b = (θ[9] / 2) / sf.beta(λ_q + 1, 5 + 1)
    
    
    # Integrate
    I_u = A_u * sf.beta(λ_u+1, K_u+1)

    I_d = A_d * sf.beta(λ_d+1, K_d+1)

    I_g = A_g1 * sf.beta(λ_g1 + 1, K_g + 1) + A_g2 * sf.beta(λ_g2 + 1, 5 + 1)

    I_q = 2 * (A_ubar + A_dbar + A_s + A_c + A_b) * sf.beta(λ_q + 1, 5 + 1)

    I_tot = I_u + I_d + I_g + I_q
    
    I_tot
end

"""
    int_xtotx(λ_u, λ_d, λ_g1, λ_g2, K_g, λ_q, θ)

Total integrated momentum density. Should equal 1.
"""
function int_xtotx(hyper_params::NamedTuple)

    hp = hyper_params

    result = int_xtotx(hp.λ_u, hp.λ_d, hp.λ_g1, hp.λ_g2,
                       hp.K_g, hp.λ_q, hp.θ)
    
    return result
end

"""
    plot_input_pdfs(hyper_params, xmin, xmax, nx)

Plot the input PDFs defined by hyper_params over 
the given x range.   
"""
function plot_input_pdfs(hyper_params::NamedTuple, xmin::Float64=1.0e-2,
                         xmax::Float64=1.0, nx::Integer=1000)

    x_grid = range(xmin, stop=xmax, length=nx)
    hp = hyper_params

    p = plot(x_grid, [x_uv_x(x, hp.λ_u, hp.θ[1]) for x in x_grid], label="x uv(x)", lw=3)
    p = plot!(x_grid, [x_dv_x(x, hp.λ_d, hp.θ[2]) for x in x_grid], label="x dv(x)", lw=3)
    p = plot!(x_grid, [x_g_x(x, hp.λ_g1, hp.λ_g2, hp.K_g, hp.θ[3], hp.θ[4]) 
                   for x in x_grid], label="x g(x)", lw=3)
    
    p = plot!(x_grid, [x_q_x(x, hp.λ_q, hp.θ[5]) for x in x_grid], label="x ubar(x)", lw=3)
    p = plot!(x_grid, [x_q_x(x, hp.λ_q, hp.θ[6]) for x in x_grid], label="x dbar(x)", lw=3)
    p = plot!(x_grid, [x_q_x(x, hp.λ_q, hp.θ[7]) for x in x_grid], label="x s(x)", lw=3)
    p = plot!(x_grid, [x_q_x(x, hp.λ_q, hp.θ[8]) for x in x_grid], label="x c(x)", lw=3)
    p = plot!(x_grid, [x_q_x(x, hp.λ_q, hp.θ[9]) for x in x_grid], label="x b(x)", lw=3)
    
    p = plot!(xlabel="x")
    p = plot!(xaxis=:log, yaxis=:log, legend=:outertopright)
    p = ylims!(1e-8, 30)

    return p
end
