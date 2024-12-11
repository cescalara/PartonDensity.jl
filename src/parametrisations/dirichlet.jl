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
@with_kw struct DirichletPDFParams{T<:Real,TV<:AbstractVector{T}} <: AbstractPDFParams
    param_type::Int = DIRICHLET_TYPE
    θ::TV
    K_u::T
    λ_u::T = (θ[1] * (K_u + 1)) / (2 - θ[1])
    K_d::T
    λ_d::T = (θ[2] * (K_d + 1)) / (1 - θ[2])
    λ_g1::T
    λ_g2::T
    K_g::T
    λ_q::T
    K_q::T
end


function xtotx(x::Real, pdf_params::DirichletPDFParams)

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

function int_xtotx_dirichlet(λ_u::Real, K_u::Real, λ_d::Real,
    K_d::Real, λ_g1::Real, λ_g2::Real, K_g::Real, λ_q::Real,
    K_q::Real, θ::AbstractArray{<:Real})

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


function get_input_pdf_func(pdf_params::DirichletPDFParams)::Function
    return function _input_pdfs(i::Integer, x::T)::T where {T<:Real}
        i = i[]
        x = x[]

        f::T = 0.0

        if i == 0
            f = parton_pdf_func(gluon, pdf_params)(x)
        elseif i == 1
            f = parton_pdf_func(u_valence, pdf_params)(x)
        elseif i == 2
            f = parton_pdf_func(d_valence, pdf_params)(x)
        elseif i == 3
            f = parton_pdf_func(u_bar, pdf_params)(x)
        elseif i == 4
            f = parton_pdf_func(d_bar, pdf_params)(x)
        elseif i == 5
            f = parton_pdf_func(s_quark, pdf_params)(x)
        elseif i == 6
            f = parton_pdf_func(s_bar, pdf_params)(x)
        elseif i == 7
            f = parton_pdf_func(c_quark, pdf_params)(x)
        elseif i == 8
            f = parton_pdf_func(c_bar, pdf_params)(x)
        elseif i == 9
            f = parton_pdf_func(b_quark, pdf_params)(x)
        elseif i == 10
            f = parton_pdf_func(b_bar, pdf_params)(x)
        end

        return f
    end
end


function parton_pdf_func end
export parton_pdf_func

function parton_pdf_func(::typeof(gluon), p::DirichletPDFParams)
    let λ_g1 = p.λ_g1, λ_g2 = p.λ_g2, K_g = p.K_g, K_q = p.K_q, w1 = p.θ[3], w2 = p.θ[4]
        return _gluon_pdf(x::Real) = x_g_x(x, λ_g1, λ_g2, K_g, K_q, w1, w2)
    end
end

function parton_pdf_func(::typeof(u_valence), p::DirichletPDFParams)
    let λ_u = p.λ_u, K_u = p.K_u
        return _u_valuence_pdf(x::Real) = x_uv_x(x, λ_u, K_u)
    end
end

function parton_pdf_func(::typeof(d_valence), p::DirichletPDFParams)
    let λ_d = p.λ_d, K_d = p.K_d
        return _d_valence_pdf(x::Real) = x_dv_x(x, λ_d, K_d)
    end
end

function parton_pdf_func(::typeof(u_bar), p::DirichletPDFParams)
    let λ_q = p.λ_q, K_q = p.K_q, θ_5 = p.θ[5]
        return _u_bar_pdf(x) = x_q_x(x, λ_q, K_q, θ_5)
    end
end

function parton_pdf_func(::typeof(d_bar), p::DirichletPDFParams)
    let λ_q = p.λ_q, K_q = p.K_q, θ_6 = p.θ[6]
        return _d_bar_pdf(x) = x_q_x(x, λ_q, K_q, θ_6)
    end
end

function parton_pdf_func(::Union{typeof(s_quark),typeof(s_bar)}, p::DirichletPDFParams)
    let λ_q = p.λ_q, K_q = p.K_q, θ_7 = p.θ[7]
        return _s_bar_pdf(x) = x_q_x(x, λ_q, K_q, θ_7)
    end
end

function parton_pdf_func(::Union{typeof(c_quark),typeof(c_bar)}, p::DirichletPDFParams)
    let λ_q = p.λ_q, K_q = p.K_q, θ_8 = p.θ[8]
        return _c_bar_pdf(x) = x_q_x(x, λ_q, K_q, θ_8)
    end
end

function parton_pdf_func(::Union{typeof(b_quark),typeof(b_bar)}, p::DirichletPDFParams)
    let λ_q = p.λ_q, K_q = p.K_q, θ_9 = p.θ[9]
        return _d_bar_pdf(x) = x_q_x(x, λ_q, K_q, θ_9)
    end
end
