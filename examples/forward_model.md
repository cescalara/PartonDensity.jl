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

# Forward model

Here, we go through an example of simulating the full forward model, from the prior definition to the expected number of events in different bins of the detector response.

```julia
using QCDNUM, PartonDensity 
using Distributions, Plots, Random, Printf, NaNMath, Parameters
pd = PartonDensity;
```

##  Define input PDFs

```julia
seed = 5
Random.seed!(seed);
dirichlet = Dirichlet([4., 4., 50., 0.5, 5., 5., 3., 2., 1.])
input_random_dirichlet = rand(dirichlet)

hyper_params = pd.PDFParameters(λ_u=0.5, λ_d=0.6, λ_g1=-0.37, λ_g2=-0.7, 
    K_g=6.0, λ_q=0.5, θ=input_random_dirichlet);
```

Plot the input PDFs

```julia
pd.plot_input_pdfs(hyper_params)
```

```julia
# Sanity check that sum = 1
pd.int_xtotx(hp.λ_u, hp.λ_d, hp.λ_g1, hp.λ_g2, hp.K_g, hp.λ_q, hp.θ) ≈ 1
```

## Define QCDNUM grids, weights and settings

```julia
QCDNUM.qcinit(-6, " ")
QCDNUM.setord(2); # 1 <=> LO, 2<=> NLO in pQCD
QCDNUM.setalf(0.118, 100.0); # α_S = 0.118, μ_R^2 = 100.0

# grid params
iosp = 3; # spline order
n_x = 100;
n_q = 50;

# x grid
xmin = Float64.([1.e-3, 0.75e0]);
iwt = Int32.([1]);
QCDNUM.gxmake(xmin, iwt, 1, n_x, iosp);
>>>>>>> main

g = qcdnum_params.grid
QCDNUM.gxmake(g.x_bounds, g.x_weights, g.x_num_bounds, g.nx, g.spline_interp);
QCDNUM.gqmake(g.qq_bounds, g.qq_weights, g.qq_num_bounds, g.nq);

# copy locally
qcdnum_x_grid = QCDNUM.gxcopy(g.nx);
qcdnum_qq_grid = QCDNUM.gqcopy(g.nq);

qc = qcdnum_params
QCDNUM.setcbt(qc.n_fixed_flav, qc.iqc, qc.iqb, qc.iqt); # 5 flavours in FFNS
iq0 = QCDNUM.iqfrmq(qcdnum_params.q0); # Get index of μ_F^2 = 100.0 = μ_R^2

nw = QCDNUM.fillwt(qcdnum_params.weight_type)
nw = QCDNUM.zmfillw()
```

## Evolve the PDFs using QCDNUM


Define input PDF function
* See https://www.nikhef.nl/~h24/qcdnum-files/doc/qcdnum170115.pdf under `evolfg`

```julia
my_func = pd.get_input_pdf_func(hyper_params)
input_pdfs = @cfunction(my_func, Float64, (Ref{Int32}, Ref{Float64}))
```

Define mapping between your input function and quark species.

```julia
map = pd.input_pdf_map;
```

```julia
eps = QCDNUM.evolfg(1, input_pdfs, map, iq0)
```

## Define necessary splines for cross section calculation

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
function _my_fun_xsec_i(ipx, ipq, first)::Float64
    
    ix = ipx[]
    iq = ipq[]
    
    xsec = _fun_xsec_i(ix, iq)
    
    return xsec
end

my_fun_xsec_i = @cfunction(_my_fun_xsec_i, Float64, (Ref{Int32}, Ref{Int32}, Ref{UInt8}))
```

```julia
_fun_xsec_i(100, 1)
```

```julia
g.nx
```

```julia
# plot
xsec_on_grid = zeros(g.nx, g.nq);

for ix = 1:g.nx
    for iq = 1:g.nq
        xsec_on_grid[ix, iq] = _fun_xsec_i(ix, iq)
    end
end

p1 = heatmap(qcdnum_x_grid, qcdnum_qq_grid, NaNMath.log10.(xsec_on_grid[:, :]'))
plot(p1, xlabel="x", ylabel="q2", 
    xaxis=:log, yaxis=:log)
```

```julia
plot(qcdnum_x_grid, NaNMath.log10.(xsec_on_grid[:, 4]), 
    label="Q2=141 (input scale)", lw=3)
plot!(qcdnum_x_grid, NaNMath.log10.(xsec_on_grid[:, 22]), label="Q2=1152", lw=3)
plot!(qcdnum_x_grid, NaNMath.log10.(xsec_on_grid[:, 35]), label="Q2=5233", lw=3)
plot!(qcdnum_x_grid, NaNMath.log10.(xsec_on_grid[:, 41]), label="Q2=10523", lw=3)
plot!(xaxis=:log, legend=:bottomleft, xlabel="x", 
    ylabel="log10(cross section spline input)", ylims=(-7, 5))
```

```julia
istx_2 = 1;
istq_2 = 2;
rscut = 370.0;
rs = 318.0;

set_lepcharge(1)
iaF_7 = QCDNUM.isp_s2make(istx_2, istq_2);
QCDNUM.ssp_uwrite(7, Float64(iaF_7));
QCDNUM.ssp_s2fill(iaF_7, my_fun_xsec_i, rscut);

set_lepcharge(-1)
iaF_8 = QCDNUM.isp_s2make(istx_2, istq_2);
QCDNUM.ssp_uwrite(8, Float64(iaF_8));
QCDNUM.ssp_s2fill(iaF_8, my_fun_xsec_i, rscut);
```

```julia
# plot spline
spline = zeros(n_x, n_q);

for ix = 1:n_x
    for iq = 1:n_q
        spline[ix, iq] = QCDNUM.dsp_funs2(iaF_7, qcdnum_x_grid[ix], 
            qcdnum_qq_grid[iq], 1)
    end
end

p1 = heatmap(qcdnum_x_grid, qcdnum_qq_grid, NaNMath.log10.(spline[:, :]'))
plot(p1, xlabel="x", ylabel="q2", 
    xaxis=:log, yaxis=:log)
```

## Integrate over the cross section spline and find expected events numbers

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

```julia
# 1 for e-p and 0 for e+p
ePp = 0;
eMp = 1; 
     
TM_eP = get_TM_elements(ePp);
TM_eM = get_TM_elements(eMp);     

K_eP = get_K_elements(ePp);
K_eM = get_K_elements(eMp);

nbins_out = size(TM_eP)[2];
```

```julia
xsec_pred_eP = zeros(nbins_out);
xsec_pred_eM = zeros(nbins_out);

for j in 1:nbins_out
    
    for i in 1:nbins
    
        xsec_pred_eP[j] += TM_eP[i, j] * (1.0/K_eP[i]) * IntXsec_eP[i];
        xsec_pred_eM[j] += TM_eM[i, j] * (1.0/K_eM[i]) * IntXsec_eM[i];
    
    end
    
end
```

```julia
xsec_pred_eM
```

```julia

```
