using QCDNUM

@with_kw struct QCDNUMGrid
    x_bounds::Vector{Float64}
    x_weights::Vector{Int32}
    x_num_bounds::Integer = size(x_bounds)[1]
    nx::Integer
    qq_bounds::Vector{Float64}
    qq_weights::Vector{Float64}
    qq_num_bounds::Integer = size(qq_bounds)[1]
    nq::Integer
    spline_interp::Integer
end

@with_kw struct QCDNUMParameters
    order::Integer
    α_S::Float64
    q0::Float64
    grid::QCDNUMGrid
    n_fixed_flav::Integer
    iqc::Integer
    iqb::Integer
    iqt::Integer
    weight_type::Integer
end
