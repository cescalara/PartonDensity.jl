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
n_points = 429 # Number of bins, should be configurable

"""
    f2_lo(x, q2)

Calculate the f2_lo structure function. To be run
after the evolution of PDFs with QCDNUM.
"""
function f2_lo(x::Float64, q2::Float64)::Float64

    # For weights
    pz = q2 / ((ZMass*ZMass + q2) * (4*Sin2ThetaW * (1 - Sin2ThetaW)))

    Au = 4.0/9.0 -2*pz*(2.0/3.0)*vu*ve + pz^2*(ve^2 + ae^2)*(vu^2 + au^2)

    Ad = 1.0/9.0 -2*pz*(-1.0/3.0)*vd*ve + pz^2*(ve^2 + ae^2)*(vd^2 + ad^2)

    # As in HERAPDF (top set to 0)
    weights = Float64.([0., Ad, Au, Ad, Au, Ad, 0., Ad, Au, Ad, Au, Ad, 0.])

    # Structure function calculation
    sum = QCDNUM.zmstfun(2, weights, x, q2, 1, 0)
    
    sum[1]
end


