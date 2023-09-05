"""
Calculation of differential cross sections
from structure functions.

Based on the Rxsecnc.h C++ code provided 
by R. Aggarwal.
"""

using QCDNUM

export QuarkCoefficients

"""
    QuarkCoefficients

Quark coefficients for structure function 
calculations.
"""
@with_kw struct QuarkCoefficients
    proup::Vector{Float64} = [0., 0., 1., 0., 1., 0., 0., 0., 1., 0., 1., 0., 0.]
    prodn::Vector{Float64} = [0., 1., 0., 1., 0., 1., 0., 1., 0., 1., 0., 1., 0.]
    valup::Vector{Float64} = [0., 0.,-1., 0.,-1., 0., 0., 0., 1., 0., 1., 0., 0.]
    valdn::Vector{Float64} = [0.,-1., 0.,-1., 0.,-1., 0., 1., 0., 1., 0., 1., 0.]
end

include("MetaData.jl")

# These globals will later be set by user
const ZMass = 91.1876
const WMass = 80.398 
const AlphaEM = 7.297352570e-03
const GFermi = 1.16637e-05 
const TopMass = 171.2 
const BottomMass = 4.20

const Vub = 41.2e-3
const Vcb = 3.93e-3

const Sin2ThetaW = 0.23127 
const Sin2ThetaC = 0.05 
const vu = 0.19164
const vd = -0.34582
const ve = -0.03746

const au = 0.5
const ad = -0.5
const ae = -0.5


export _fun_xsec_i, get_input_xsec_func

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
    rxsecnc_xq2_i(charge, md, x, q2, F2, xF3, FL)

Reduced cross section for single x, q2.
"""
function rxsecnc_xq2_i(charge::Int, md::MetaData, x::Float64, q2::Float64, F2::Float64, xF3::Float64, FL::Float64)::Float64

    rxsec = -1.0 #?
    y = 0.04  #?
    
    y = q2 / md.sqrtS / md.sqrtS / x
    Y_plus = 1 + (1 - y)^2
    Y_minus = 1 - (1 - y)^2

    if (charge == 1)
        
        rxsec =  F2 - (Y_minus / Y_plus)*xF3 - (y^2 / Y_plus) * FL
        
    elseif (charge == -1)

        rxsec =  F2 + (Y_minus / Y_plus)*xF3 - (y^2 / Y_plus) * FL
        
    end

    rxsec
end

"""
    rxsecnc_xq2(x, q2, md)

Reduced cross section for all bins.
"""
function rxsecnc_xq2(chage::Int, md::MetaData, x_bin_cen::Array{Float64}, q2_bin_cen::Array{Float64})::Array{Float64}

    n_bins = length(x_bin_cen)

    # Check same number of x and q2 vals
    if (n_bins != length(q2_bin_cen))

        throw(DimensionMismatch("x_bin_cen and q2_bin_cen must have same length."))

    end

    # Calculate for each bin
    rxsec = zeros(n_bins)

    for i = 1:n_bins
        
        rxsec[i] = rxsecnc_xq2_i(charge, md, x_bin_cen[i], q2_bin_cen[i])
        
    end
    
    rxsec
end

"""
    nc_propagator(md, q2, x)
"""
function nc_propagator(md::MetaData, q2::Float64, x::Float64)::Float64

    y = q2 / md.sqrtS / md.sqrtS / x
    Yplus = 1 + (1 - y)^2
    alpha = AlphaEM
    conversion_factor = 0.3894e9;  # GeV^2 -> pb^2
    
    conversion_factor*2*π*Yplus*alpha^2 / (x * q2^2)
end

"""
    dd_xsecnc_xq2_i(charge, md,  x, q2)

Double differential cross section for single 
x and q2 values. 
NB: modifications needed to include pol and order.
"""
function dd_xsecnc_xq2_i(charge::Int, md::MetaData, x::Float64, q2::Float64, F2::Float64, xF3::Float64, FL::Float64)::Float64

    dd_xsec = -1.0
    
    dd_xsec = nc_propagator(md, q2, x) * rxsecnc_xq2_i(charge, md, x, q2, F2, xF3, FL)
    
    return dd_xsec
end

"""
    dd_xsecnc_xq2(charge, x_bin_cen, q2_bin_cen)

Double differential cross section for all x and 
q2 bins.
NB: modifications needed to include pol and order.
"""
function dd_xsecnc_xq2(charge::Int, md::MetaData, x_bin_cen::Array{Float64}, q2_bin_cen::Array{Float64})::Array{Float64}

    n_bins = length(x_bin_cen)

    # Check same number of x and q2 vals
    if (n_bins != length(q2_bin_cen))

        throw(DimensionMismatch("x_bin_cen and q2_bin_cen must have same length."))

    end

    xsec = zeros(n_bins)
    
    for i = 1:n_bins
        
        xsec[i] = dd_xsecnc_xq2_i(charge, md, x_bin_cen[i], q2_bin_cen[i])
        
    end
    
    xsec
end

"""
    _fun_xsec_i(ix iq)

Input function for cross section spline.
Must be wrapped for interface to SPLINT.
"""
function _fun_xsec_i(charge::Int,md::MetaData, ix, iq)::Float64
    
    # get q2 and x values
    q2 = QCDNUM.qfrmiq(iq);
    x = QCDNUM.xfrmix(ix);

    # get spline addresses
    iF2up = Int32(QCDNUM.dsp_uread(1));
    iF2dn = Int32(QCDNUM.dsp_uread(2));
    iF3up = Int32(QCDNUM.dsp_uread(3));
    iF3dn = Int32(QCDNUM.dsp_uread(4));
    iFLup = Int32(QCDNUM.dsp_uread(5));
    iFLdn = Int32(QCDNUM.dsp_uread(6));
    
    # structure function calculation
    pz = q2 / ((ZMass*ZMass+q2) * (4*(Sin2ThetaW) * (1-Sin2ThetaW)));

    Au = 4.0/9.0 -2*pz*(2.0/3.0)*(vu)*(ve) + pz*pz*(ve*ve+ae*ae)*(vu*vu+au*au);
    Ad = 1.0/9.0 -2*pz*(-1.0/3.0)*(vd)*(ve) + pz*pz*(ve*ve+ae*ae)*(vd*vd+ad*ad);
    Bu = -2*(2.0/3.0)*au*ae*pz + 4*au*ae*vu*ve*pz*pz;
    Bd = -2*(-1.0/3.0)*ad*ae*pz + 4*ad*ae*vd*ve*pz*pz;
    
    F2 = Au * QCDNUM.dsp_funs2(iF2up, x, q2, 1) + Ad * QCDNUM.dsp_funs2(iF2dn, x, q2, 1);
    xF3 = Bu * QCDNUM.dsp_funs2(iF3up, x, q2, 1) + Bd * QCDNUM.dsp_funs2(iF3dn, x, q2, 1);
    FL = Au * QCDNUM.dsp_funs2(iFLup, x, q2, 1) + Ad * QCDNUM.dsp_funs2(iFLdn, x, q2, 1);
    
    xsec = dd_xsecnc_xq2_i(charge, md, x, q2, F2, xF3, FL);
    
    return  xsec;
end

function get_input_xsec_func(charge::Int, md::MetaData)

if ( charge == 1 )
    funcp = function _my_fun_xsec_ip(ipx, ipq, first)::Float64
    
        ix = ipx[]
        iq = ipq[]
        mycharge::Int = 1
        
        xsec = _fun_xsec_i(mycharge, md, ix, iq)
    
        return xsec
    end
        return funcp
end
if (charge == -1 ) 

    funcm = function _my_fun_xsec_im(ipx, ipq, first)::Float64
    
        ix = ipx[]
        iq = ipq[]
        mycharge::Int = -1
        
        xsec = _fun_xsec_i(mycharge, md, ix, iq)
    
        return xsec
    end

    return funcm
end

end
