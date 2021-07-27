---
jupyter:
  jupytext:
    text_representation:
      extension: .md
      format_name: markdown
      format_version: '1.3'
      jupytext_version: 1.11.0
  kernelspec:
    display_name: Julia 1.6.1
    language: julia
    name: julia-1.6
---

```julia
using QCDNUM, PartonDensity 
using Distributions, Plots, Random, Printf, NaNMath
pd = PartonDensity;
```

```julia
iord = 3
nfin = 0
as0 = 0.364
r20 = 2.0

def = Float64.([
        0., 0., 0., 0., 0.,-1., 0., 1., 0., 0., 0., 0., 0.,
        0., 0., 0., 0., 0., 1., 0., 0., 0., 0., 0., 0., 0.,    
        0., 0., 0., 0.,-1., 0., 0., 0., 1., 0., 0., 0., 0.,   
        0., 0., 0., 0., 1., 0., 0., 0., 0., 0., 0., 0., 0.,   
        0., 0., 0.,-1., 0., 0., 0., 0., 0., 1., 0., 0., 0.,    
        0., 0., 0., 1., 0., 0., 0., 0., 0., 0., 0., 0., 0.,    
        0., 0.,-1., 0., 0., 0., 0., 0., 0., 0., 1., 0., 0.,    
        0., 0., 1., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,    
        0.,-1., 0., 0., 0., 0., 0., 0., 0., 0., 0., 1., 0.,    
        0., 1., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
        0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
        0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
]);

xmin = Float64.([1e-3, 0.2, 0.4, 0.6, 0.75])
iwt = Int32.([1, 2, 4, 8, 16]);

ngx = 5
nxin = 100
iosp = 3
qq = Float64.([1.0, 2.0, 25, 3e4])
wt = Float64.([1.0, 1.0, 1.0 , 1.0])
ngq = 4 
nqin = 50
```

```julia
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

What do the test PDFs look like?

```julia
x_grid = range(0.01, stop=1, length=500)

plot(x_grid, [func(0, x) for x in x_grid], lw=3, label="x glu(x)")
plot!(x_grid, [func(1, x) for x in x_grid], lw=3, label="x dval(x)")
plot!(x_grid, [func(2, x) for x in x_grid], lw=3, label="x dbar(x)")
plot!(x_grid, [func(3, x) for x in x_grid], lw=3, label="x uval(x)")
plot!(x_grid, [func(4, x) for x in x_grid], lw=3, label="x ubar(x)")
plot!(x_grid, [func(5, x) for x in x_grid], lw=3, label="x sval(x)")
plot!(x_grid, [func(6, x) for x in x_grid], lw=3, label="x sbar(x)")
plot!(xaxis=:log, xlabel="x")
```

Evolve the test PDF over the defined qq grid.

```julia
QCDNUM.qcinit("/usr/local/lib/libQCDNUM.dylib", -6, " ")
nx = QCDNUM.gxmake(xmin, iwt, ngx, nxin, iosp)
nq = QCDNUM.gqmake(qq, wt, ngq, nqin)
nw = QCDNUM.fillwt(1)
nw = QCDNUM.zmfillw()
QCDNUM.setord(iord)
QCDNUM.setalf(as0, r20)
iqc = QCDNUM.iqfrmq(qq[2])
iqb = QCDNUM.iqfrmq(qq[3])
QCDNUM.setcbt(0, iqc, iqb, 999)
iq0 = QCDNUM.iqfrmq(qq[1])
iq1 = nq
qq1 = QCDNUM.qfrmiq(iq1)
eps = QCDNUM.evolfg(1, func_c, def, iq0);
```

Plot slices in F2

```julia
# Copy grids locally
qcdnum_x_grid = QCDNUM.gxcopy(nx);
qcdnum_qq_grid = QCDNUM.gqcopy(nq);
```

```julia
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

function test_f2_lo(x::Float64, q2::Float64)::Float64

    # For weights
    pz = q2 / ((ZMass^2 + q2) * (4*Sin2ThetaW * (1 - Sin2ThetaW)))

    Au = 4.0/9.0 #-2*pz*(2.0/3.0)*vu*ve + pz^2*(ve^2 + ae^2)*(vu^2 + au^2)

    Ad = 1.0/9.0 #-2*pz*(-1.0/3.0)*vd*ve + pz^2*(ve^2 + ae^2)*(vd^2 + ad^2)

    # As in HERAPDF (top set to 0)
    weights = [0., Ad, Au, Ad, Au, Ad, 0., Ad, Au, Ad, Au, Ad, 0.]

    # Structure function calculation
    output = QCDNUM.zmstfun(2, weights, [x], [q2], 1, 0)
    
    output[1]
end
```

```julia
Nx = size(qcdnum_x_grid)[1]
Nq = size(qcdnum_qq_grid)[1]
F_test = zeros(Nx, Nq)

for ix = 1:Nx
    for iq = 1:Nq
        F_test[ix, iq] = test_f2_lo(qcdnum_x_grid[ix], qcdnum_qq_grid[iq]) 
    end
end
```

```julia
qcdnum_qq_grid[49]
```

```julia
# Slices in Q2
plot(qcdnum_x_grid, F_test[:, 1], label="Q2=1 (input scale)", lw=3)
#plot!(qcdnum_x_grid, F_test[:, 25], label="Q2=160", lw=3)
#plot!(qcdnum_x_grid, F_test[:,40], label="Q2=3,700", lw=3)
#plot!(qcdnum_x_grid, F_test[:,49], label="Q2=24,000", lw=3)
plot!(xaxis=:log, legend=:topright, xlabel="x")
```

```julia

```
