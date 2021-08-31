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
    display_name: Julia 1.6.1
    language: julia
    name: julia-1.6
---

```julia
using QCDNUM, PartonDensity 
using Distributions, Plots, Random, Printf, NaNMath
pd = PartonDensity;
```

### Select hyperparameters

```julia
seed = 5
#seed = 100;
#seed = 200; 
#seed = 300;
#seed = 42
#seed = rand(10:10000)
```

```julia
#Random.seed!(seed);
#dirichlet = Dirichlet([5., 4., 1, 0.5, 0.3, 0.2, 0.1, 0.1, 0.1])
#input_random_dirichlet = rand(dirichlet)

#hp = NamedTuple{(:λ_u, :λ_d, :λ_g1, :λ_g2, :K_g, :λ_q, :θ)}((0.99, 0.99, 0.7, -0.7, 1.1, -0.5, input_random_dirichlet));
```

```julia
Random.seed!(seed);
dirichlet = Dirichlet([4., 4., 50., 0.5, 5., 5., 3., 2., 1.])
input_random_dirichlet = rand(dirichlet)

hp = NamedTuple{(:λ_u, :λ_d, :λ_g1, :λ_g2, :K_g, :λ_q, :θ)}((0.5, 0.6, -0.37, -0.7, 6., -0.5, input_random_dirichlet));
```

```julia
input_random_dirichlet
```

### For simple comparison

```julia
using SpecialFunctions, Roots
const sf = SpecialFunctions
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

### Plot input PDFs

```julia
x_grid = range(0.001, stop=1, length=1000)

plot(x_grid, [pd.x_uv_x(x, hp.λ_u, hp.θ[1]) for x in x_grid], label="x uv(x)", lw=3)
plot!(x_grid, [pd.x_dv_x(x, hp.λ_d, hp.θ[2]) for x in x_grid], label="x dv(x)", lw=3)
plot!(x_grid, [pd.x_g_x(x, hp.λ_g1, hp.λ_g2, hp.K_g, hp.θ[3], hp.θ[4]) 
        for x in x_grid], label="x g(x)", lw=3)
plot!(x_grid, [pd.x_q_x(x, hp.λ_q, hp.θ[5]) for x in x_grid], label="x ubar(x)", lw=3)
plot!(x_grid, [pd.x_q_x(x, hp.λ_q, hp.θ[6]) for x in x_grid], label="x dbar(x)", lw=3)
plot!(x_grid, [pd.x_q_x(x, hp.λ_q, hp.θ[7]) for x in x_grid], label="x s(x)", lw=3)
plot!(x_grid, [pd.x_q_x(x, hp.λ_q, hp.θ[8]) for x in x_grid], label="x c(x)", lw=3)
plot!(x_grid, [pd.x_q_x(x, hp.λ_q, hp.θ[9]) for x in x_grid], label="x b(x)", lw=3)
plot!(xlabel="x")
ylims!(1e-8, 30)
plot!(xaxis=:log, yaxis=:log, legend=:outertopright)
#plot!(xaxis=:log, legend=:outertopright)
```

```julia
# Sanity check that sum = 1
pd.int_xtotx(hp.λ_u, hp.λ_d, hp.λ_g1, hp.λ_g2, hp.K_g, hp.λ_q, hp.θ)
```

### Evolve using QCDNUM


Make x-qq grid
* Can some of this be hidden from the user?
* Can we silence the output? - not so important as this will be called once per fit.

```julia
xmin = Float64.([1e-3, 0.2, 0.4, 0.6, 0.75]) # Where grid density changes
iwt = Int32.([1, 2, 4, 8, 16]) # Weights for more grid points at higher x
nxin = 200 # Request 100 x grid points
iord = 3 # Quadratic interpolation

qlim = Float64.([100, 3e4]) # Limits of qq grid
wt = Float64.([1e0, 1e0]) # Weights for even grid points
nqin = 50 # Request 50 qq grid points
```

```julia
QCDNUM.qcinit("/usr/local/lib/libQCDNUM.dylib", -6, " ")
nx = QCDNUM.gxmake(xmin, iwt, size(xmin,1), nxin, iord)
nq = QCDNUM.gqmake(qlim, wt, size(qlim,1), nqin)

itype = 1 # Unpolarised
nw = QCDNUM.fillwt(itype)
nw = QCDNUM.zmfillw()

@printf("Generated grid with size nx = %i, nq = %i.\n", nx, nq)
@printf("Filled weight tables.\n")
```

```julia
# Copy grids locally
qcdnum_x_grid = QCDNUM.gxcopy(nx);
qcdnum_qq_grid = QCDNUM.gqcopy(nq);
```

Set input parameters 
* Could some of this be hidden from the user?

```julia
QCDNUM.setord(2); # 1 <=> LO, 2<=> NLO in pQCD
QCDNUM.setalf(0.118, 100.0); # α_S = 0.118, μ_R^2 = 100.0
QCDNUM.setcbt(5, 1, 1, 1); # 5 flavours in FFNS
iq0 = QCDNUM.iqfrmq(100.0); # Get index of μ_F^2 = 100.0 = μ_R^2
```

Pass input PDF function
* See https://www.nikhef.nl/~h24/qcdnum-files/doc/qcdnum170115.pdf under `evolfg`

```julia code_folding=[]
function _input_pdfs(i, x)::Float64
    i = i[]
    x = x[]
    
    f = 0.0
    
    # gluon
    if (i == 0)
        f = pd.x_g_x(x, hp.λ_g1, hp.λ_g2, hp.K_g, hp.θ[3], hp.θ[4]) 
    end
  
    # u valence
    if (i == 1)
        f = pd.x_uv_x(x, hp.λ_u, hp.θ[1])
    end
    
    # d valence
    if (i == 2)
        f = pd.x_dv_x(x, hp.λ_d, hp.θ[2])
    end
    
    # ubar
    if (i == 3)
        f = pd.x_q_x(x, hp.λ_q, hp.θ[5])
    end
    
    # dbar
    if (i == 4)
        f = pd.x_q_x(x, hp.λ_q, hp.θ[6])
    end
    
    # s and sbar
    if (i == 5) || (i == 6)
        f = pd.x_q_x(x, hp.λ_q, hp.θ[7])
    end
    
    # c and cbar
    if (i == 7) || (i == 8)
        f = pd.x_q_x(x, hp.λ_q, hp.θ[8])
    end
    
    # d and dbar
    if (i == 9) || (i == 10)
        f = pd.x_q_x(x, hp.λ_q, hp.θ[9])
    end
    
    return f
end

input_pdfs = @cfunction(_input_pdfs, Float64, (Ref{Int32}, Ref{Float64}))
```

This acts as a mapping between your input function and quark species.

```julia
#               tb  bb  cb  sb  ub  db  g   d   u   s   c   b   t
map = Float64.([0., 0., 0., 0.,-1., 0., 0., 0., 1., 0., 0., 0., 0., # 1 # U valence
                0., 0., 0., 0., 0.,-1., 0., 1., 0., 0., 0., 0., 0., # 2 # D valence
                0., 0., 0., 0., 1., 0., 0., 0., 0., 0., 0., 0., 0., # 3 # u sea
                0., 0., 0., 0., 0., 1., 0., 0., 0., 0., 0., 0., 0., # 4 # d sea
                0., 0., 0., 0., 0., 0., 0., 0., 0., 1., 0., 0., 0., # 5 # s
                0., 0., 0., 1., 0., 0., 0., 0., 0., 0., 0., 0., 0., # 6 # sbar
                0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 1., 0., 0., # 7 # c
                0., 0., 1., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., # 8 # cbar
                0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 1., 0., # 9 # b
                0., 1., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., # 10 # bbar
                0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., # 11
                0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.]); # 12
```

```julia
iset1 = 1                                       
jtype = 10*iset1+itype
eps = QCDNUM.evolfg(jtype, input_pdfs, map, iq0)
#eps = QCDNUM.evolfg(jtype, func_c, def, iq0) - need to evolve to correct scales...
```

Make function to pass structure function evaluations to SPLINT

```julia
function _func_to_integrate(ix, iq, first)::Float64
    
    # Deref pointers
    ix = ix[] 
    iq = iq[]

    # Get x, qq values
    x = QCDNUM.xfrmix(ix)
    q2 = QCDNUM.qfrmiq(iq)
    
    # Get double differential cross section    
    return dd_xsecnc_xq2_i(x, q2)    
end

func_to_integrate = @cfunction(_func_to_integrate, Float64, (Ref{Int32}, Ref{Int32}, Ref{UInt8}))
```

```julia
Nx = size(qcdnum_x_grid)[1]
Nq = size(qcdnum_qq_grid)[1]
F = zeros(Nx, Nq);
F_test = zeros(Nx, Nq)

for ix = 1:Nx
    for iq = 1:Nq
        F[ix, iq] = _func_to_integrate(ix, iq, false)
        F_test[ix, iq] = f2_lo(qcdnum_x_grid[ix], qcdnum_qq_grid[iq]) 
    end
end
```

### Input vs interpolated PDFs

```julia
# Find interpolated PDFs with sumfxq
ichk = 1

wd = Float64.([0., 0., 0., 0., 0., 0., 0., 1., 0., 0., 0., 0., 0.]);
wu = Float64.([0., 0., 0., 0.,0., 0., 0., 0., 1., 0., 0., 0., 0.]);
wdbar = Float64.([0., 0., 0., 0., 0.,1., 0., 0., 0., 0., 0., 0., 0.]);
wubar = Float64.([0., 0., 0., 0.,1., 0., 0., 0., 0., 0., 0., 0., 0.]);
ws = Float64.([0., 0., 0., 0., 0.,0., 0., 0., 0., 1., 0., 0., 0.]);
wc = Float64.([0., 0., 0., 0.,0., 0., 0., 0., 0., 0., 1., 0., 0.]);
wsbar = Float64.([0., 0., 0., 1., 0.,0., 0., 0., 0., 0., 0., 0., 0.]);
wcbar = Float64.([0., 0., 1., 0.,0., 0., 0., 0., 0., 0., 0., 0., 0.]);
wb = Float64.([0., 0., 0., 0.,0., 0., 0., 0., 0., 0., 0., 1., 0.]);
wbbar = Float64.([0., 1., 0., 0., 0.,0., 0., 0., 0., 0., 0., 0., 0.]);

q_ = 100.0
x_grid = range(0.01, stop=1, length=500)

u = [QCDNUM.sumfxq(iset1,wu,1,x,q_,ichk) for x in x_grid];
ubar = [QCDNUM.sumfxq(iset1,wubar,1,x,q_,ichk) for x in x_grid];
d = [QCDNUM.sumfxq(iset1,wd,1,x,q_,ichk) for x in x_grid];
dbar = [QCDNUM.sumfxq(iset1,wdbar,1,x,q_,ichk) for x in x_grid];
s = [QCDNUM.sumfxq(iset1,ws,1,x,q_,ichk) for x in x_grid];
sbar = [QCDNUM.sumfxq(iset1,wsbar,1,x,q_,ichk) for x in x_grid];
c = [QCDNUM.sumfxq(iset1,wc,1,x,q_,ichk) for x in x_grid];
cbar = [QCDNUM.sumfxq(iset1,wcbar,1,x,q_,ichk) for x in x_grid];
b = [QCDNUM.sumfxq(iset1,wb,1,x,q_,ichk) for x in x_grid]
bbar = [QCDNUM.sumfxq(iset1,wbbar,1,x,q_,ichk) for x in x_grid];
g = [QCDNUM.sumfxq(iset1, wu, 0, x, q_, ichk) for x in x_grid];
```

```julia
# Plot comparison
#plot(x_grid, [pd.x_uv_x(x, hp.λ_u, hp.θ[1]) for x in x_grid], label="x u(x) input", 
#    lw=3, alpha=1)
#plot(x_grid, [pd.x_dv_x(x, hp.λ_d, hp.θ[2]) for x in x_grid], label="x d(x)", lw=3)
plot(x_grid, [pd.x_g_x(x, hp.λ_g1, hp.λ_g2, hp.K_g, hp.θ[3], hp.θ[4]) 
        for x in x_grid], label="x g(x)", lw=3)
#plot(x_grid, [pd.x_q_x(x, hp.λ_q, hp.θ[5]) for x in x_grid], label="x ubar(x)", lw=3)
#plot(x_grid, [pd.x_q_x(x, hp.λ_q, hp.θ[6]) for x in x_grid], label="x dbar(x)", lw=3)
#plot(x_grid, [pd.x_q_x(x, hp.λ_q, hp.θ[7]) for x in x_grid], label="x s(x)", lw=3)
#plot(x_grid, [pd.x_q_x(x, hp.λ_q, hp.θ[8]) for x in x_grid], label="x c(x)", lw=3)
#plot(x_grid, [pd.x_q_x(x, hp.λ_q, hp.θ[9]) for x in x_grid], label="x b(x)", lw=3)

plot!(x_grid, g, label="x g(x) interp", lw=3, linestyle=:dash)
plot!(xlabel="x")
#ylims!(1e-8, 10)
plot!(xaxis=:log, legend=:outertopright)
```

### F2 and differential cross-section

x bounds = (1e-3, 1)
q2 bounds = (1e2, 3e4)

```julia
# What does this function look like?
sel = qcdnum_qq_grid .> 300; # For useful bins
p1 = heatmap(qcdnum_x_grid, qcdnum_qq_grid, F_test[:, :]')
plot(p1, xlabel="x", ylabel="q2", 
    xaxis=:log, yaxis=:log)
```

```julia
qcdnum_qq_grid[41]
```

```julia
# Slices in Q2
plot(qcdnum_x_grid, F_test[:, 1], label="Q2=100 (input scale)")
plot!(qcdnum_x_grid, F_test[:, 22], label="Q2=1026")
plot!(qcdnum_x_grid, F_test[:,35], label="Q2=5234")
plot!(qcdnum_x_grid, F_test[:,41], label="Q2=10523")
plot!(xaxis=:log, legend=:topleft)
```

```julia
# Differential cross section
p1 = heatmap(qcdnum_x_grid, qcdnum_qq_grid, NaNMath.log10.(F[:, :]'))
plot(p1, xlabel="x", ylabel="q2", title="log10(Diff. cross section)", 
    xaxis=:log, yaxis=:log)
```

Make spline object

```julia
QCDNUM.ssp_spinit(1000)
iasp = QCDNUM.isp_s2make(5, 5)

# Fill the spline and set no kinematic limit
QCDNUM.ssp_s2fill(iasp, func_to_integrate, 0.0)

# Set user nodes to avoid numerical issues
#x_nodes = Vector(range(3e-3, stop=0.9, length=20))
#qq_nodes = Vector(range(3e2, stop=2e3, length=10))
#nx = length(x_nodes)
#nq = length(qq_nodes)
#iasp = QCDNUM.isp_s2user(x_nodes, nx, qq_nodes, nq)

# Print spline summary
QCDNUM.ssp_nprint(iasp)
```

Define binning

```julia
Nx = size(qcdnum_x_grid)[1]
Nq = size(qcdnum_qq_grid)[1]
spline = zeros(Nx, Nq);

for ix = 1:Nx
    for iq = 1:Nq
        spline[ix, iq] = QCDNUM.dsp_funs2(iasp, qcdnum_x_grid[ix], 
            qcdnum_qq_grid[iq], 1)
    end
end
```

Visualise the spline

```julia
p1 = heatmap(qcdnum_x_grid, qcdnum_qq_grid, NaNMath.log10.(spline[:, :]'))
plot(p1, xlabel="x", ylabel="q2", title="log10(Spline)", 
    xaxis=:log, yaxis=:log)
```

Integrated cross section

```julia
integ_xsec = zeros(length(qcdnum_x_grid)-1, length(qcdnum_qq_grid)-1);

for ix = 1:length(qcdnum_x_grid)-1
    for iq = 1:length(qcdnum_qq_grid)-1
        integ_xsec[ix, iq] = QCDNUM.dsp_ints2(iasp, qcdnum_x_grid[ix],
            qcdnum_x_grid[ix+1], 
            qcdnum_qq_grid[iq], qcdnum_qq_grid[iq+1])
    end
end
```

```julia
p1 = heatmap(qcdnum_x_grid[1:end-1], qcdnum_qq_grid[1:end-1], 
    NaNMath.log10.(integ_xsec[:, :]'))
plot(p1, xlabel="x", ylabel="q2", title="log10(Int. cross section)", 
    xaxis=:log, yaxis=:log)
```

```julia

```
