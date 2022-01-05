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

```julia
using QCDNUM, PartonDensity 
using Distributions, Plots, Random, Printf, NaNMath
pd = PartonDensity;
```

## Define PDFs

```julia
seed = 5
Random.seed!(seed);
dirichlet = Dirichlet([4., 4., 50., 0.5, 5., 5., 3., 2., 1.])
input_random_dirichlet = rand(dirichlet)

hyper_params = pd.PDFParameters(λ_u=0.5, λ_d=0.6, λ_g1=-0.37, λ_g2=-0.7, 
    K_g=6.0, λ_q=0.5, θ=input_random_dirichlet);
```

```julia
pd.plot_input_pdfs(hyper_params)
```

```julia
# Sanity check that sum = 1
pd.int_xtotx(hyper_params) ≈ 1
```

## Define grid, weights and evolve PDFs

```julia
# grid parameters
iosp = 3; # spline order
n_x = 100;
n_q = 50;

QCDNUM.qcinit(-6, " ")
QCDNUM.setord(2); # 1 <=> LO, 2<=> NLO in pQCD
QCDNUM.setalf(0.118, 100.0); # α_S = 0.118, μ_R^2 = 100.0

# x grid
xmin = Float64.([1.e-3, 0.2e0, 0.4e0, 0.6e0, 0.75e0]);
iwt = Int32.([1, 2, 4, 8, 16]);
ngx = 5
QCDNUM.gxmake(xmin, iwt, 1, n_x, iosp);

# mu2 grid
qarr = Float64.([1.e2, 3.e4]);
warr = Float64.([1., 1.]);
ngq = 2;
QCDNUM.gqmake(qarr, warr, ngq, n_q);

# copy locally
qcdnum_x_grid = QCDNUM.gxcopy(n_x);
qcdnum_qq_grid = QCDNUM.gqcopy(n_q);

QCDNUM.setcbt(5, 1, 1, 1); # 5 flavours in FFNS
iq0 = QCDNUM.iqfrmq(100.0); # Get index of μ_F^2 = 100.0 = μ_R^2

itype = 1 # Unpolarised
nw = QCDNUM.fillwt(itype)
nw = QCDNUM.zmfillw()
```

```julia
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
# For splines
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
```

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

```julia

```
