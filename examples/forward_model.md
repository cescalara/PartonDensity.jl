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
using Distributions, Plots, Random, Printf
pd = PartonDensity
```

### Select hyperparameters

```julia
Random.seed!(5);
dirichlet = Dirichlet([5., 4., 1, 0.5, 0.3, 0.2, 0.1, 0.1, 0.1])

hp = NamedTuple{(:λ_u, :λ_d, :λ_g1, :λ_g2, :K_g, :λ_q, :θ)}((0.7, 0.5, 0.7, -0.7, 1.0, 
        -0.5, rand(dirichlet)))
```

### Plot input PDFs

```julia
x_grid = range(0, stop=1, length=50)

plot(x_grid, [pd.x_uv_x(x, hp.λ_u, hp.θ[1]) for x in x_grid], label="x uv(x)", lw=3)
plot!(x_grid, [pd.x_dv_x(x, hp.λ_d, hp.θ[2]) for x in x_grid], label="x dv(x)", lw=3)
plot!(x_grid, [pd.x_g_x(x, hp.λ_g1, hp.λ_g2, hp.K_g, hp.θ[3], hp.θ[4]) 
        for x in x_grid], label="x g2(x)", lw=3)
plot!(x_grid, [pd.x_q_x(x, hp.λ_q, hp.θ[5]) for x in x_grid], label="x ubar(x)", lw=3)
plot!(x_grid, [pd.x_q_x(x, hp.λ_q, hp.θ[6]) for x in x_grid], label="x dbar(x)", lw=3)
plot!(x_grid, [pd.x_q_x(x, hp.λ_q, hp.θ[7]) for x in x_grid], label="x s(x)", lw=3)
plot!(x_grid, [pd.x_q_x(x, hp.λ_q, hp.θ[8]) for x in x_grid], label="x c(x)", lw=3)
plot!(x_grid, [pd.x_q_x(x, hp.λ_q, hp.θ[9]) for x in x_grid], label="x b(x)", lw=3)
plot!(xlabel="x")
```

### Evolve using QCDNUM


Make x-qq grid
* Can some of this be hidden from the user?
* Can we silence the output?

```julia
xmin = Float64.([1e-5, 0.2, 0.4, 0.6, 0.75]) # Where grid density changes
iwt = Int32.([1, 2, 4, 8, 16]) # Weights for more grid points at higher x
nxin = 100 # Request 100 x grid points
iord = 3 # Quadratic interpolation

qlim = Float64.([2e0, 1e4]) # Limits of qq grid
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
QCDNUM.setord(2); # NLO in pQCD
QCDNUM.setalf(0.364, 2.0); # α_S = 0.364, μ_R^2 = 2.0
QCDNUM.setcbt(6, 1, 1, 1); # 6 flavours in FFNS
iq0 = QCDNUM.iqfrmq(2.0); # Get index of μ_F^2 = 2.0 = μ_R^2
```

Pass input PDF function
* How best to split this up in `i` and `j` as explained in the QCDNUM manual?
* See https://www.nikhef.nl/~h24/qcdnum-files/doc/qcdnum170115.pdf under `evolfg`

```julia
function func(i, x)::Float64
    i = i[]
    x = x[]
    
    f = 0.0
    
    # gluon
    if (i == 0)
        f = pd.x_g_x(x, hp.λ_g1, hp.λ_g2, hp.K_g, hp.θ[3], hp.θ[4]) 
    end
  
    # uv ?
    if (i == 1)
        f = pd.x_uv_x(x, hp.λ_u, θ[1])
    end
    
    # ud ?
    if (i == 2)
        f = pd.x_ud_x(x, hp.λ_d, θ[2])
    end
    
    # etc...
    
end
```

```julia
# This acts as a mapping between your input function and quark species
# From examples
def = Float64.([0., 0., 0., 0., 0.,-1., 0., 1., 0., 0., 0., 0., 0., 
                0., 0., 0., 0.,-1., 0., 0., 0., 1., 0., 0., 0., 0., 
                0., 0., 0.,-1., 0., 0., 0., 0., 0., 1., 0., 0., 0.,
                0., 0., 0., 0., 0., 1., 0., 0., 0., 0., 0., 0., 0., 
                0., 0., 0., 0., 1., 0., 0., 0., 0., 0., 0., 0., 0., 
                0., 0., 0., 1., 0., 0., 0., 0., 0., 0., 0., 0., 0., 
                0., 0.,-1., 0., 0., 0., 0., 0., 0., 0., 1., 0., 0., 
                0., 0., 1., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 
                0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 
                0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 
                0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 
                0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.]);
```
