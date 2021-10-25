---
jupyter:
  jupytext:
    formats: ipynb,md
    text_representation:
      extension: .md
      format_name: markdown
      format_version: '1.3'
      jupytext_version: 1.11.0
  kernelspec:
    display_name: Julia 1.7.0-rc2
    language: julia
    name: julia-1.7
---

## Test new forward model implementation

Use a simple test PDF from the QCDNUM example jobs to test and benchmark the forward model. Based on `testzmstf.cc` by R. Aggarwal.

```julia
using QCDNUM, PartonDensity 
using Distributions, Plots, Random, Printf, NaNMath
pd = PartonDensity;
```

Define mapping for input PDF in QCDNUM, as well as structure of outputs and splines.

```julia code_folding=[31]
               # tb  bb  cb  sb  ub  db   g   d   u   s   c   b   t
               # -6  -5  -4  -3  -2  -1   0   1   2   3   4   5   6
pdfin = Float64.([0., 0., 0., 0., 0.,-1., 0., 1., 0., 0., 0., 0., 0.,   # dval
                  0., 0., 0., 0.,-1., 0., 0., 0., 1., 0., 0., 0., 0.,   # uval
                  0., 0., 0.,-1., 0., 0., 0., 0., 0., 1., 0., 0., 0.,   # sval
                  0., 0., 0., 0., 0., 1., 0., 0., 0., 0., 0., 0., 0.,   # dbar
                  0., 0., 0., 0., 1., 0., 0., 0., 0., 0., 0., 0., 0.,   # ubar
                  0., 0., 0., 1., 0., 0., 0., 0., 0., 0., 0., 0., 0.,   # sbar
                  0., 0.,-1., 0., 0., 0., 0., 0., 0., 0., 1., 0., 0.,   # cval
                  0., 0., 1., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,   # cbar
                  0.,-1., 0., 0., 0., 0., 0., 0., 0., 0., 0., 1., 0.,   # bval
                  0., 1., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,   # bbar
                  -1., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 1.,  # tval
                  1., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.]); # tbar
    
             # tb  bb  cb  sb  ub  db   g   d   u   s   c   b   t
             # -6  -5  -4  -3  -2  -1   0   1   2   3   4   5   6
dnv = Float64.([0., 0., 0., 0., 0.,-1., 0., 1., 0., 0., 0., 0., 0.]);
upv = Float64.([0., 0., 0., 0.,-1., 0., 0., 0., 1., 0., 0., 0., 0.]);
del = Float64.([0., 0., 0., 0.,-1., 1., 0., 0., 0., 0., 0., 0., 0.]);
uds = Float64.([0., 0., 0., 2., 2., 2., 0., 0., 0., 0., 0., 0., 0.]);
pro = Float64.([4., 1., 4., 1., 4., 1., 0., 1., 4., 1., 4., 1., 4.]) / 9;
duv = Float64.([0., 0., 0., 0.,-1.,-1., 0., 1., 1., 0., 0., 0., 0.]);
val = Float64.([-1.,-1.,-1.,-1.,-1.,-1., 0., 1., 1., 1., 1., 1., 1.]);

proup = Float64.([0., 0., 1., 0., 1., 0., 0., 0., 1., 0., 1., 0., 0.]);
prodn = Float64.([0., 1., 0., 1., 0., 1., 0., 1., 0., 1., 0., 1., 0.]);
valup = Float64.([0.,0.,-1.,0.,-1.,0., 0., 0., 1., 0., 1., 0., 0.]);
valdn = Float64.([0.,-1.,0.,-1.,0.,-1., 0., 1., 0., 1., 0., 1., 0.]);

# input pdf
function func(ipdf, x)::Float64
    i = ipdf[]
    xb = x[]
    adbar = 0.1939875
    f = 0
    if (i == 0) 
        ag = 1.7
        f = ag * xb^-0.1 * (1.0-xb)^5.0
    end
    if (i == 1)
        ad = 3.064320
        f = ad * xb^0.8 * (1.0-xb)^4.0
    end
    if (i == 2)
        au = 5.107200
        f = au * xb^0.8 * (1.0-xb)^3.0
    end
    if (i == 3) 
        f = 0.0
    end
    if (i == 4)
        f = adbar * xb^-0.1 * (1.0-xb)^6.0
    end
    if (i == 5) 
        f = adbar * xb^-0.1 * (1.0-xb)^6.0 * (1.0-xb)
    end
    if (i == 6)
        xdbar = adbar * xb^-0.1 * (1.0-xb)^6.0
        xubar = adbar * xb^-0.1 * (1.0-xb)^6.0 * (1.0-xb)
        f = 0.2 * (xdbar + xubar)
    end
    if (i == 7)
        f = 0.0
    end
    if (i == 8)
        f = 0.0
    end
    if (i == 9)
        f = 0.0
    end
    if (i == 10)
        f = 0.0
    end
    if (i == 11)
        f = 0.0
    end
    if (i == 12)
        f = 0.0
    end
    return f
end

func_c = @cfunction(func, Float64, (Ref{Int32}, Ref{Float64}))
```

Define more key inputs.

```julia
# cbt quark masses
hqmass = Float64.([1.43, 5., 80.]);

# LO/NLO/NNLO
iord = 2;

# VFNS
nfix = 0;

# define alpha_s
alf = 0.35e0;
q2a = 2.0e0;

# more grid parameters
iosp = 3; # spline order
n_x = 100;
n_q = 50;
q0 = 2.e0;
```

Initialise QCDNUM and make grids and wieght tables. Then evolve the PDF.

```julia
QCDNUM.qcinit("/usr/local/lib/libQCDNUM.dylib", -6, " ")
QCDNUM.setord(iord)
QCDNUM.setalf(alf, q2a)

# x grid
xmin = Float64.([1.e-3, 0.2e0, 0.4e0, 0.6e0, 0.75e0]);
iwt = Int32.([1, 2, 4, 8, 16]);
ngx = 5
QCDNUM.gxmake(xmin, iwt, 1, n_x, iosp);

# mu2 grid
qarr = Float64.([1.e0, 2.e0, 25.e0, 1.e5]);
warr = Float64.([4., 4., 4., 4.]);
ngq = 4;
QCDNUM.gqmake(qarr, warr, ngq, n_q);

# VFNS
iq0 = QCDNUM.iqfrmq(qarr[1]);
iqc = QCDNUM.iqfrmq(hqmass[1] * hqmass[1]);
iqb = QCDNUM.iqfrmq(hqmass[2] * hqmass[2]);
iqt = QCDNUM.iqfrmq(hqmass[3] * hqmass[3]);
QCDNUM.setcbt(nfix, iqc, iqb, iqt);

# weights
itype = 1; # unpolarised
nw = QCDNUM.fillwt(itype);
nw = QCDNUM.zmfillw();

# evolution
QCDNUM.evolfg(1, func_c, pdfin, iq0);
```

Define structure function splines

```julia
nuser = 10
QCDNUM.ssp_spinit(nuser);
istx = 5;
istq = 10;
ia = QCDNUM.isp_s2make(istx, istq);

xnd = QCDNUM.ssp_unodes(ia, 100, 0);
qnd = QCDNUM.ssp_vnodes(ia, 100, 0);
QCDNUM.ssp_nprint(ia);
QCDNUM.ssp_erase(ia);

#  edit nodes
xnd[91] = 0.13e0;
xnd[92] = 0.16e0;
xnd[93] = 0.33e0;
```

```julia
# set nodes and fill spline with structure function
iaFLup = QCDNUM.isp_s2user(xnd, 100, qnd, 100);
QCDNUM.ssp_s2f123(iaFLup, 1, proup, 1, 0.0);

iaF2up = QCDNUM.isp_s2user(xnd, 100, qnd, 100);
QCDNUM.ssp_s2f123(iaF2up, 1, proup, 2, 0.0);

iaF3up = QCDNUM.isp_s2user(xnd, 100, qnd, 100);
QCDNUM.ssp_s2f123(iaF3up, 1, valup, 3, 0.0);
    
iaFLdn = QCDNUM.isp_s2user(xnd, 100, qnd, 100);
QCDNUM.ssp_s2f123(iaFLdn, 1, prodn, 1, 0.0);

iaF2dn = QCDNUM.isp_s2user(xnd, 100, qnd, 100);
QCDNUM.ssp_s2f123(iaF2dn, 1, prodn, 2, 0.0);

iaF3dn = QCDNUM.isp_s2user(xnd, 100, qnd, 100);
QCDNUM.ssp_s2f123(iaF3dn, 1, valdn, 3, 0.0);

# store spline addresses
QCDNUM.ssp_uwrite(1, Float64(iaF2up));
QCDNUM.ssp_uwrite(2, Float64(iaF2dn));
QCDNUM.ssp_uwrite(3, Float64(iaF3up));
QCDNUM.ssp_uwrite(4, Float64(iaF3dn));
QCDNUM.ssp_uwrite(5, Float64(iaFLup));
QCDNUM.ssp_uwrite(6, Float64(iaFLdn));
```

Define cross section and cross section splines.


```julia code_folding=[25, 48, 58, 68]
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

function my_rxsecnc_xq2_i(x::Float64, q2::Float64, F2::Float64, 
        xF3::Float64, FL::Float64)::Float64

    rxsec = -1.0 #?
    y = 0.04  #?
    
    y = q2 / sqrt_s / sqrt_s / x
    Y_plus = 1 + (1 - y)^2
    Y_minus = 1 - (1 - y)^2

    if (Lepcharge == 1)
        
        rxsec =  F2 - (Y_minus / Y_plus)*xF3 - (y^2 / Y_plus) * FL
        
    elseif (Lepcharge == -1)

        rxsec =  F2 + (Y_minus / Y_plus)*xF3 - (y^2 / Y_plus) * FL
        
    end

    rxsec
end

function my_nc_propagator(q2::Float64, x::Float64)::Float64

    y = q2 / sqrt_s / sqrt_s / x
    Yplus = 1 + (1 - y)^2
    alpha = AlphaEM
    conversion_factor = 0.3894e9;  # GeV^2 -> pb^2
    
    conversion_factor*2*π*Yplus*alpha^2 / (x * q2^2)
end

function my_dd_xsecnc_xq2_i(x::Float64, q2::Float64, F2::Float64, 
        xF3::Float64, FL::Float64)::Float64

    dd_xsec = -1.0
    
    dd_xsec = my_nc_propagator(q2, x) * my_rxsecnc_xq2_i(x, q2, F2, xF3, FL)
    
    return dd_xsec
end

function _my_fun_xsec_i(ipx, ipq, first)::Float64
    
    ix = ipx[]
    iq = ipq[]
    
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
    F3 = Bu * QCDNUM.dsp_funs2(iF3up, x, q2, 1) + Bd * QCDNUM.dsp_funs2(iF3dn, x, q2, 1);
    FL = Au * QCDNUM.dsp_funs2(iFLup, x, q2, 1) + Ad * QCDNUM.dsp_funs2(iFLdn, x, q2, 1);
    
    xsec = my_dd_xsecnc_xq2_i(x, q2, F2, F3, FL);

    return xsec;    
end

my_fun_xsec_i = @cfunction(_my_fun_xsec_i, Float64, (Ref{Int32}, Ref{Int32}, Ref{UInt8}))
```

```julia
# cross section spline
istx_2 = 1;
istq_2 = 2;
rscut = 370.0;
rs = 318.0;

Lepcharge = 1;
iaF_7 = QCDNUM.isp_s2make(istx_2, istq_2);
QCDNUM.ssp_uwrite(7, Float64(iaF_7));
QCDNUM.ssp_s2fill(iaF_7, my_fun_xsec_i, rscut);

Lepcharge = -1;
iaF_8 = QCDNUM.isp_s2make(istx_2, istq_2);
QCDNUM.ssp_uwrite(8, Float64(iaF_8));
QCDNUM.ssp_s2fill(iaF_8, my_fun_xsec_i, rscut);
```

Calculate the integrated cross sections.

```julia
nbins = size(xbins_M_begin)[1]
IntXsec_eP = zeros(nbins);
IntXsec_eM = zeros(nbins);
for i in 1:nbins
    IntXsec_eP[i] = QCDNUM.dsp_ints2(iaF_7, xbins_M_begin[i], xbins_M_end[i], 
        q2bins_M_begin[i], q2bins_M_end[i], 318., 4);
    IntXsec_eM[i] = QCDNUM.dsp_ints2(iaF_8, xbins_M_begin[i], xbins_M_end[i], 
        q2bins_M_begin[i],q2bins_M_end[i], 318., 4);
end     
```

Calculate expected event number in bins

```julia
# 1 for e-p and 0 for e+p
ePp = 0;
eMp = 1; 
     
TM_eP = get_TM_elements(ePp);
TM_eM = get_TM_elements(eMp);     

K_eP = get_K_elements(ePp);
K_eM = get_K_elements(eMp);

nbins_out = size(TM_eP)[2];
xsec_pred_eP = zeros(nbins_out);
xsec_pred_eM = zeros(nbins_out);

for j in 1:nbins_out
    
    for i in 1:nbins
    
        xsec_pred_eP[j] += TM_eP[i, j] * (1.0/K_eP[i]) * IntXsec_eP[i];
        xsec_pred_eM[j] += TM_eM[i, j] * (1.0/K_eM[i]) * IntXsec_eM[i];
    
    end
    
    xsec_pred_eP[j] *= get_L_data(ePp);
    xsec_pred_eM[j] *= get_L_data(eMp);
    
end
```

```julia

```
