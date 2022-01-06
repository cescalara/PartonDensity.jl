using QCDNUM

export QCDNUMGrid, QCDNUMParameters
export SPLINTParameters

"""
    QCDNUMGrid

Struct for holding the QCDNUM grid parameters.
"""
@with_kw struct QCDNUMGrid
    "x boundaries of grid"
    x_min::Vector{Float64}
    "x grid weights"
    x_weights::Vector{Int32}
    "number of x grid boundaries"
    x_num_bounds::Integer = size(x_min)[1]
    "number of x grid points"
    nx::Integer
    "qq boundaries of grid"
    qq_bounds::Vector{Float64}
    "qq grid weights"
    qq_weights::Vector{Float64}
    "number of qq grid boundaries"
    qq_num_bounds::Integer = size(qq_bounds)[1]
    "number of qq grid points"
    nq::Integer
    "degree of spline interpolation used"
    spline_interp::Integer = 3
end

"""
   QCDNUMParameters

Struct for holding all QCDNUM Parameters. 
"""
@with_kw struct QCDNUMParameters
    "order of evolution in pQCD"
    order::Integer
    "coupling constant at starting scale"
    α_S::Float64
    "starting scale"
    q0::Float64
    "QCDNUMGrid"
    grid::QCDNUMGrid
    "number of fixed flavours in FFNS"
    n_fixed_flav::Integer = 5
    "charm threshold as index of qq grid (VFNS only)"
    iqc::Integer = 1
    "bottom threshold as index of qq grid (VFNS only)"
    iqb::Integer = 1
    "top threshold as index of qq grid (VFNS only)"
    iqt::Integer = 1
    "weight table type (1<=>unpolarised)"
    weight_type::Integer = 1
end

"""
    SplineAddresses

Lookup table for addresses of different 
structure function splines.
"""
@with_kw struct SplineAddresses
    F2up::Integer = 1
    F2dn::Integer = 2
    F3up::Integer = 3
    F3dn::Integer = 4
    FLup::Integer = 5
    FLdn::Integer = 6
    F_eP::Integer = 7
    F_eM::Integer = 8
end

"""
    SPLINTParameters

Struct for storage of parameters used
with SPLINT package of QCDNUM.
"""
@with_kw struct SPLINTParameters
    "number of words in memory for user space"
    nuser::Integer = 10
    "number of steps in x"
    nsteps_x::Integer = 5
    "number of steps in qq"
    nsteps_q::Integer = 10
    "number of nodes in x"
    nnodes_x::Integer = 100
    "number of nodes in qq"
    nnodes_q::Integer = 100
    "rs constraint"
    rs::Float64 = 318.0
    "cut on rs"
    rscut::Float64 = 370.0
    "spline addresses"
    spline_addresses::SplineAddresses = SplineAddresses()
end
