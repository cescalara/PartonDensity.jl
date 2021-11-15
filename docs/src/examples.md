# Examples

## Forward model

Here, we go through an example of simulating the full forward model, from the prior definition to the expected number of events in different bins of the detector response.

```@example 1
using QCDNUM, PartonDensity 
using Distributions, Plots, Random, Printf, NaNMath
pd = PartonDensity;
```

### Defining the input PDFS

We start by defining the input PDFs based on our prior model.

```@example 1
seed = 5
Random.seed!(seed);
dirichlet = Dirichlet([4., 4., 50., 0.5, 5., 5., 3., 2., 1.])
input_random_dirichlet = rand(dirichlet)

hp = NamedTuple{(:λ_u, :λ_d, :λ_g1, :λ_g2, :K_g, :λ_q, :θ)}((0.5, 0.6, -0.37, -0.7, 6., -0.5, input_random_dirichlet));
```

We can plot these to understand what they look like.

```@example 1
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
```

The intgeral over these components should equal 1 - let's check that.

```@example 1
pd.int_xtotx(hp.λ_u, hp.λ_d, hp.λ_g1, hp.λ_g2, hp.K_g, hp.λ_q, hp.θ) ≈ 1
```

### Defining QCDNUM grids, weights and settings

Next, we use the QCDNUM.jl interface to set things up.

```@example 1
QCDNUM.qcinit(-6, " ")
QCDNUM.setord(2); # 1 <=> LO, 2<=> NLO in pQCD
QCDNUM.setalf(0.118, 100.0); # α_S = 0.118, μ_R^2 = 100.0

# grid params
iosp = 3; # spline order
n_x = 100;
n_q = 50;

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

@printf("Generated grid with size nx = %i, nq = %i.\n", n_x, n_q)
@printf("Filled weight tables.\n")
```

### Evolve the PDFs using QCDNUM

Once this is done, we can define a function to pass out input PDFs to QCDNUM and perform the evolution.

```@example 1
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

# we need to wrap this to pass to Fortran
input_pdfs = @cfunction(_input_pdfs, Float64, (Ref{Int32}, Ref{Float64}))

# mapping between your input function and quark species
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


```@example 1
eps = QCDNUM.evolfg(1, input_pdfs, map, iq0)
```


### Defining the necessary splines for the cross section calculation

We first make splines of the structure functions for fast cross-section calculation.

```@example 1
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

```@example 1
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

Similarly to the input PDFs, we also have to define a function to calculate the cross section.

```@example 1
function _my_fun_xsec_i(ipx, ipq, first)::Float64
    
    ix = ipx[]
    iq = ipq[]
    
    xsec = _fun_xsec_i(ix, iq)
    
    return xsec
end

# wrap to pass to Fortran
my_fun_xsec_i = @cfunction(_my_fun_xsec_i, Float64, (Ref{Int32}, Ref{Int32}, Ref{UInt8}))

# fill spline
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

### Integrate over the cross section spline and find expected event numbers

```@example 1
nbins = size(xbins_M_begin)[1]
IntXsec_eP = zeros(nbins);
IntXsec_eM = zeros(nbins);
for i in 1:nbins
    IntXsec_eP[i] = QCDNUM.dsp_ints2(iaF_7, xbins_M_begin[i], xbins_M_end[i], 
        q2bins_M_begin[i], q2bins_M_end[i], 318., 4);
    IntXsec_eM[i] = QCDNUM.dsp_ints2(iaF_8, xbins_M_begin[i], xbins_M_end[i], 
        q2bins_M_begin[i],q2bins_M_end[i], 318., 4);
end 

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
    
end
```

```@example 1
xsec_pred_eM
```
