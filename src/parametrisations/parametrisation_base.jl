import SpecialFunctions as SF

export xtotx, int_xtotx
export plot_input_pdfs, get_input_pdf_func
export input_pdf_map

"""
    abstract type PartonDensity.AbstractPDFParams

Abstract type for any PDF parametrisation.
"""
abstract type AbstractPDFParams end

"""
    x_uv_x(x, λ_u, K_u)

Momentum density of u valence component.
Beta function 
    A_u x^λ_u (1-x)^K_u
A_u is set by λ_u and K_u.
"""
function x_uv_x(x::Real, λ_u::Real, K_u::Real)

    A_u = 2 / SF.beta(λ_u, K_u + 1)

    return A_u * x^λ_u * (1 - x)^K_u

end

"""
    x_dv_x(x, λ_d, K_d)

Momentum density of d valence component.
Beta function
    A_d x^λ_d (1 - x)^K_d
A_d is set by λ_d and K_d.
"""
function x_dv_x(x::Real, λ_d::Real, K_d::Real)

    A_d = 1 / SF.beta(λ_d, K_d + 1)

    return A_d * x^λ_d * (1 - x)^K_d

end

"""
        x_g_x(x, λ_g1, λ_g2, K_g, K_q, w1, w2)

Momentum density of gluon component.
    A_g1 x^λ_g1 (1 - x)^K_g + A_g2 x^λ_g2 (1 - x)^K_q
Amplitudes are set by weights `w1` and `w2`.
"""
function x_g_x(x::Real, λ_g1::Real, λ_g2::Real, K_g::Real,
    K_q::Real, w1::Real, w2::Real)

    A_g1 = w1 / SF.beta(λ_g1 + 1, K_g + 1)

    x_g1_x = A_g1 * x^λ_g1 * (1 - x)^K_g

    A_g2 = w2 / SF.beta(λ_g2 + 1, K_q + 1)

    x_g2_x = A_g2 * x^λ_g2 * (1 - x)^K_q

    return x_g1_x + x_g2_x

end


"""
    x_q_x(x, λ_q, K_q, w)

Momentum density of non-valence
quark component.
    A_q x^λ_q (1 - x)^K_q 
Amplitude is set by weight `w`
"""
function x_q_x(x::Real, λ_q::Real, K_q::Real, w::Real)

    A_q = (w / 2) / SF.beta(λ_q + 1, K_q + 1)

    return A_q * x^λ_q * (1 - x)^K_q

end

"""
    x_total_x(x, pdf_params)

Total momentum density.
"""
function xtotx end

"""
    int_xtotx(pdf_params)

Total integrated momentum density. Should equal 1.
"""
function int_xtotx end

function plot_input_pdfs end
"""
    plot_input_pdfs(pdf_params; xmin, xmax, nx)

Plot the input PDFs defined by hyper_params over 
the given x range.   
"""

function get_input_pdf_func end

"""
    input_pdf_map

The relevant mapping for use with QCDNUM and `get_input_pdf_func`.
"""
const input_pdf_map = Float64.([
#   tb  bb  cb  sb  ub  db   g   d   u   s   c   b   t
     0,  0,  0,  0, -1 , 0,  0,  0,  1,  0,  0,  0,  0, # 1 # U valence
     0,  0,  0,  0,  0, -1,  0,  1,  0,  0,  0,  0,  0, # 2 # D valence
     0,  0,  0,  0,  1,  0,  0,  0,  0,  0,  0,  0,  0, # 3 # u sea
     0,  0,  0,  0,  0,  1,  0,  0,  0,  0,  0,  0,  0, # 4 # d sea
     0,  0,  0,  0,  0,  0,  0,  0,  0,  1,  0,  0,  0, # 5 # s
     0,  0,  0,  1,  0,  0,  0,  0,  0,  0,  0,  0,  0, # 6 # sbar
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  1,  0,  0, # 7 # c
     0,  0,  1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, # 8 # cbar
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  1,  0, # 9 # b
     0,  1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, # 10 # bbar
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, # 11
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0  # 12
])
