"""
Calculation of differential cross sections
from structure functions.

Based on the Rxsecnc.h C++ code provided 
by R. Aggarwal.
"""

using QCDNUM

# These globals will later be set by user
ZMass = 91.1876
WMass = 80.398 
AlphaEM = 7.297352570e-03
GFermi = 1.16637e-05 
TopMass = 171.2 
BottomMass = 4.20

Vub = 41.2e-3
Vcb = 3.93e-3

Sin2ThetaW = 0.23127 
Sin2ThetaC = 0.05 
vu = 0.19164
vd = -0.34582
ve = -0.03746

au = 0.5
ad = -0.5
ae = -0.5

# Variables for Rxsecnc_xq2 
sqrt_s = 318.1 # Should be configurable for calculation of y
Lepcharge = 1 # Should be configurable for calculation of Y

export dd_xsecnc_xq2_i

"""
    f2_lo(x, q2)

Calculate the f2_lo structure function term. 
To be run after the evolution of PDFs with QCDNUM.
"""
function f2_lo(x::Float64, q2::Float64)::Float64

    # For weights
    pz = q2 / ((ZMass^2 + q2) * (4*Sin2ThetaW * (1 - Sin2ThetaW)))

    Au = 4.0/9.0 -2*pz*(2.0/3.0)*vu*ve + pz^2*(ve^2 + ae^2)*(vu^2 + au^2)

    Ad = 1.0/9.0 -2*pz*(-1.0/3.0)*vd*ve + pz^2*(ve^2 + ae^2)*(vd^2 + ad^2)

    # As in HERAPDF (top set to 0)
    weights = [0., Ad, Au, Ad, Au, Ad, 0., Ad, Au, Ad, Au, Ad, 0.]

    # Structure function calculation
    output = QCDNUM.zmstfun(2, weights, [x], [q2], 1, 0)
    
    output[1]
end

"""
    xf3_lo(x, q2)

Calculate the xf3_lo structure function term.
To be run after the evolution of PDFs with QCDNUM.
"""
function xf3_lo(x::Float64, q2::Float64)::Float64

    # For weights
    pz = q2 / ((ZMass^2 + q2) * (4*Sin2ThetaW * (1 - Sin2ThetaW)))

    Bu = -2*(2.0/3.0)*au*ae*pz + 4*au*ae*vu*ve*pz^2
    
    Bd = -2*(-1.0/3.0)*ad*ae*pz + 4*ad*ae*vd*ve*pz^2
    
    # weights = [-Bu, -Bd, -Bu, -Bd, -Bu, -Bd, 0.,  Bd, Bu, Bd, Bu, Bd, Bu]

    # As in HERAPDF (top set to 0)
    weights = [0., -Bd, -Bu, -Bd, -Bu, -Bd, 0.,  Bd, Bu, Bd, Bu, Bd, 0.] 

    # Structure function calculation
    output = QCDNUM.zmstfun(3, weights, [x], [q2], 1, 0)

    output[1]
end

"""
    fl_lo(x, q2)

Calculate the fl_lo structure function term.
To be run after the evolution of PDFs with QCDNUM.
"""
function fl_lo(x::Float64, q2::Float64)::Float64

    # For weights
    pz = q2 / ((ZMass^2 + q2) * (4*Sin2ThetaW * (1 - Sin2ThetaW)))

    Au = 4.0/9.0 -2*pz*(2.0/3.0)*vu*ve + pz^2*(ve^2 + ae^2)*(vu^2 + au^2)
    
    Ad = 1.0/9.0 -2*pz*(-1.0/3.0)*vd*ve + pz^2(ve^2 + ae^2)*(vd^2 + ad^2)

    # As in HERAPDF (top set to 0)
    weights = [0., Ad, Au, Ad, Au, Ad, 0.,  Ad, Au, Ad, Au, Ad, 0.]

    output = QCDNUM.zmstfun(1, weights, [x], [q2], 1, 0)

    output[1]
end

"""
    rxsecnc_xq2_i(x, q2)

Reduced cross section for single x, q2.
"""
function rxsecnc_xq2_i(x::Float64, q2::Float64)::Float64

    rxsec = -1.0 #?
    y = 0.04  #?
    y = q2 / sqrt_s / sqrt_s / x
    
    Y_plus = 1 + (1 - y)^2
    Y_minus = 1 - (1 - y)^2

    F2 = f2_lo(x, q2)
    xF3 = xf3_lo(x, q2)
    FL = fl_lo(x, q2)
  
    if (Lepcharge == 1)
        
        rxsec =  F2 - (Y_minus / Y_plus)*xF3 - (y^2 / Y_plus) * FL
        
    elseif (Lepcharge == -1)

        rxsec =  F2 + (Y_minus / Y_plus)*xF3 - (y^2 / Y_plus) * FL
        
    end

    rxsec
end

"""
    rxsecnc_xq2(x, q2)

Reduced cross section for all bins.
"""
function rxsecnc_xq2(x_bin_cen::Array{Float64}, q2_bin_cen::Array{Float64})::Array{Float64}

    n_bins = length(x_bin_cen)

    # Check same number of x and q2 vals
    if (n_bins != length(q2_bin_cen))

        throw(DimensionMismatch("x_bin_cen and q2_bin_cen must have same length."))

    end

    # Calculate for each bin
    rxsec = zeros(n_bins)

    for i = 1:n_bins
        
        rxsec[i] = rxsecnc_xq2_i(x_bin_cen[i], q2_bin_cen[i])
        
    end
    
    rxsec
end

"""
    nc_propagator(q2, x)
"""
function nc_propagator(q2::Float64, x::Float64)::Float64

    y = q2 / sqrt_s / sqrt_s / x
    Yplus = 1 + (1 - y)^2
    alpha = AlphaEM
    conversion_factor = 0.3894e9;  # GeV^2 -> pb^2
    
    conversion_factor*2*π*Yplus*alpha^2 / (x * q2^2)
end

"""
    dd_xsecnc_xq2_i(x, q2)

Double differential cross section for single 
x and q2 values. 
NB: modifications needed to include pol and order.
"""
function dd_xsecnc_xq2_i(x::Float64, q2::Float64)

    dd_xsec = -1.0
    
    dd_xsec = nc_propagator(q2, x) * rxsecnc_xq2_i(x, q2)
    
    dd_xsec
end

"""
    dd_xsecnc_xq2(x_bin_cen, q2_bin_cen)

Double differential cross section for all x and 
q2 bins.
NB: modifications needed to include pol and order.
"""
function dd_xsecnc_xq2(x_bin_cen::Array{Float64},
                       q2_bin_cen::Array{Float64})::Array{Float64}

    n_bins = length(x_bin_cen)

    # Check same number of x and q2 vals
    if (n_bins != length(q2_bin_cen))

        throw(DimensionMismatch("x_bin_cen and q2_bin_cen must have same length."))

    end

    xsec = zeros(n_bins)
    
    for i = 1:n_bins
        
        xsec[i] = dd_xsecnc_xq2_i(x_bin_cen[i], q2_bin_cen[i])
        
    end
    
  xsec
end


