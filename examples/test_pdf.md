---
jupyter:
  jupytext:
    text_representation:
      extension: .md
      format_name: markdown
      format_version: '1.3'
      jupytext_version: 1.13.0
  kernelspec:
    display_name: Julia 1.7.0-beta4
    language: julia
    name: julia-1.7
---

## Test PDF

Checking the implementation of a test PDF. Based on the `batune**.f` testjobs.

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
wt = Float64.([1.0, 1.0, 1.0, 1.0])
ngq = 4 
nqin = 50
```

```julia code_folding=[]
function func(ipdf, x)::Float64
    i = ipdf[]
    xb = x[]
    
    f = 0.0
    if (i == 0) 
        ag = 1.7
        f = ag * xb^-0.1 * (1.0-xb)^5.0
    end
    if (i == 1)
        ad = 3.064320
        f = ad * xb^0.8 * (1.0-xb)^4.0
    end
    if (i == 2)
        adbar = 0.1939875
        f = adbar * xb^-0.1 * (1.0-xb)^6.0
    end
    if (i == 3) 
        au = 5.107200
        f = au * xb^0.8 * (1.0-xb)^3.0
    end
    if (i == 4)
        adbar = 0.1939875
        f = adbar * xb^-0.1 * (1.0-xb)^6.0 * (1.0-xb)
    end
    if (i == 5) 
        f = 0.0
    end
    if (i == 6)
        adbar = 0.1939875
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
x_grid = range(0.001, stop=1, length=500)
ids = [0, 1, 2, 3, 4, 5, 6]
labels = ["x glu(x)", "x dval(x)", "x dbar(x)", "x uval(x)", "x ubar(x)", 
    "x sval(x)", "x sbar(x)"]
plot()
for (id, label) in zip(ids, labels)
    plot!(x_grid, [func(id, x) for x in x_grid], lw=3, label=label)
end
plot!(xaxis=:log, xlabel="x")
```

Evolve the test PDF over the defined qq grid. We use VFNS NNLO evolution.

```julia
#Adjust the location  and name of Fortran QCDNUM library on your system, e.g. use /usr/local/lib/libQCDNUM.dylib or /usr/lib64/libQCDNUM.so.
QCDNUM.qcinit("/usr/local/lib/libQCDNUM.dylib", -6, " ")
nx = QCDNUM.gxmake(xmin, iwt, ngx, nxin, iosp)
nq = QCDNUM.gqmake(qq, wt, ngq, nqin)
nw = QCDNUM.fillwt(1)
nw = QCDNUM.zmfillw()
QCDNUM.setord(iord)
QCDNUM.setalf(as0, r20)
iqc = QCDNUM.iqfrmq(qq[2])
iqb = QCDNUM.iqfrmq(qq[3])
QCDNUM.setcbt(nfin, iqc, iqb, 999)
iq0 = QCDNUM.iqfrmq(qq[1])
iq1 = nq
qq1 = QCDNUM.qfrmiq(iq1)
eps = QCDNUM.evolfg(1, func_c, def, iq0);
```

Plot input scale vs. evolved PDFs

```julia
# Copy grids locally
qcdnum_x_grid = QCDNUM.gxcopy(nx);
qcdnum_qq_grid = QCDNUM.gqcopy(nq);

f_ids = [(0, 0, 0), (1, -1, -1), (-1, 0, 0), (2, -2, -1), 
    (-2, 0, 0), (3, -3, -1), (-3, 0, 0)];
```

```julia
# Compare with inputs
plot()
for (f_id, label) in zip(f_ids, labels)
    id1, id2, s = f_id
    pdf = ([QCDNUM.fvalij(1, id1, ix, 1, 1) for ix=1:nx] 
        .+ s*[QCDNUM.fvalij(1, id2, ix, 1, 1) for ix=1:nx])
    plot!(qcdnum_x_grid, pdf, label=label, lw=3)     
end
for (id, label) in zip(ids, labels)
    plot!(x_grid, [func(id, x) for x in x_grid], lw=3, label=label, 
        color="black", linestyle=:dash)
end
plot!(xaxis=:log, xlabel="x", title="Input scale, Q2=1")
```

```julia
iq = QCDNUM.iqfrmq(100.0)
plot()
for (f_id, label) in zip(f_ids, labels)
    id1, id2, s = f_id
    pdf = ([QCDNUM.fvalij(1, id1, ix, iq, 1) for ix=1:nx] 
        .+ s*[QCDNUM.fvalij(1, id2, ix, iq, 1) for ix=1:nx])
    plot!(qcdnum_x_grid, pdf, label=label, lw=3)     
end
plot!(xaxis=:log, legend=:topright, xlabel="x", title="Evolved scale, Q2=100")
```

```julia
plot()
for (f_id, label) in zip(f_ids, labels)
    id1, id2, s = f_id
    pdf = ([QCDNUM.fvalij(1, id1, ix, nq, 1) for ix=1:nx] 
        .+ s*[QCDNUM.fvalij(1, id2, ix, nq, 1) for ix=1:nx])
    plot!(qcdnum_x_grid, pdf, label=label, lw=3)     
end
plot!(xaxis=:log, legend=:topright, xlabel="x", title="Evolved scale, Q2=3e4")
```

Compare with `batune00.f` at Q2 = 100 GeV.

```julia
iq = QCDNUM.iqfrmq(100.0)
sel = 2:3
plot()
for (f_id, label) in zip(f_ids[sel], labels[sel])
    id1, id2, s = f_id
    pdf = ([QCDNUM.fvalij(1, id1, ix, iq, 1) for ix=1:nx] 
        .+ s*[QCDNUM.fvalij(1, id2, ix, iq, 1) for ix=1:nx])
    plot!(qcdnum_x_grid, pdf, label=label, lw=3)     
end
plot!(xaxis=:log, legend=:topright, xlabel="x", title="Evolved scale, Q2=100")
```

Plot slices in F2

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
F_test = zeros(nx, nq)

for ix = 1:nx
    for iq = 1:nq
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

## Fit parametrisation to the test pdf

Start with the evolved PDFs at Q2=100.

```julia
using SpecialFunctions, Roots
sf = SpecialFunctions
```

```julia
iq = QCDNUM.iqfrmq(100.0)
sel = 4:5
plot()
for (f_id, label) in zip(f_ids[sel], labels[sel])
    id1, id2, s = f_id
    pdf = ([QCDNUM.fvalij(1, id1, ix, iq, 1) for ix=1:nx] 
        .+ s*[QCDNUM.fvalij(1, id2, ix, iq, 1) for ix=1:nx])
    plot!(qcdnum_x_grid, pdf, label=label, lw=3)     
end
plot!(xaxis=:log, legend=:topright, xlabel="x", title="Evolved scale, Q2=100")
```

Compare with your parametrisation

```julia
seed = 5
Random.seed!(seed);
dirichlet = Dirichlet([4., 4., 50., 0.5, 5., 5., 3., 2., 1.])
input_random_dirichlet = rand(dirichlet)

hp = NamedTuple{(:λ_u, :λ_d, :λ_g1, :λ_g2, :K_g, :λ_q, :θ)}((0.5, 0.6, -0.37, -0.7, 6., -0.5, input_random_dirichlet));
```

```julia
input_random_dirichlet
```

```julia
x_grid = qcdnum_x_grid
plot()
plot!(x_grid, [pd.x_uv_x(x, hp.λ_u, hp.θ[1]) for x in x_grid], label="x uv(x)", lw=3)
#plot!(x_grid, [pd.x_dv_x(x, hp.λ_d, hp.θ[2]) for x in x_grid], label="x dv(x)", lw=3)
plot!(x_grid, [pd.x_g_x(x, hp.λ_g1, hp.λ_g2, hp.K_g, hp.θ[3], hp.θ[4]) 
        for x in x_grid], label="x g(x)", lw=3)
plot!(x_grid, [pd.x_q_x(x, hp.λ_q, hp.θ[5]) for x in x_grid], label="x ubar(x)", lw=3)
#plot!(x_grid, [pd.x_q_x(x, hp.λ_q, hp.θ[6]) for x in x_grid], label="x dbar(x)", lw=3)
#plot!(x_grid, [pd.x_q_x(x, hp.λ_q, hp.θ[7]) for x in x_grid], label="x s(x)", lw=3)
#plot!(x_grid, [pd.x_q_x(x, hp.λ_q, hp.θ[8]) for x in x_grid], label="x c(x)", lw=3)
#plot!(x_grid, [pd.x_q_x(x, hp.λ_q, hp.θ[9]) for x in x_grid], label="x b(x)", lw=3)
plot!(xlabel="x")
plot!(xaxis=:log, legend=:outertopright)
```

```julia
w = hp.θ[1]
f(K_u) = w - 2 * (sf.beta(hp.λ_u + 1, K_u + 1) / sf.beta(hp.λ_u, K_u + 1))    
K_u = find_zero(f, 1)
```

```julia
w = hp.θ[2]
f(K_d) = w - (sf.beta(hp.λ_d + 1, K_d + 1) / sf.beta(hp.λ_d, K_d + 1))    
K_d = find_zero(f, 1)
```

```julia
pd.int_xtotx(hp.λ_u, hp.λ_d, hp.λ_g1, hp.λ_g2, hp.K_g, hp.λ_q, hp.θ)
```

```julia

```
