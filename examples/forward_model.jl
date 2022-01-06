# # Forward model
# 
# Here, we go through an example of simulating the full forward model,
# from the prior definition to the expected number of events in different
# bins of the detector response.

using QCDNUM, PartonDensity 
using Plots, Printf, NaNMath, Parameters

# ## Define input PDFs

weights = [50., 0.5, 5., 5., 3., 2., 1.]
hyper_params = PDFParameters(λ_u=0.5, K_u=4.0, λ_d=0.6, K_d=6.0, λ_g1=-0.37, λ_g2=-0.7,
                             K_g=6.0, λ_q=0.5, seed=5, weights=weights);

# Plot the input PDFs

plot_input_pdfs(hyper_params)

# Sanity check that sum = 1

int_xtotx(hyper_params) ≈ 1

# ## Define QCDNUM grids, weights and settings

grid = QCDNUMGrid(x_bounds=[1.0e-3, 0.75], x_weights=[1], nx=100,
                  qq_bounds=[1.0e2, 3.0e4], qq_weights=[1.0, 1.0], nq=50, spline_interp=3)

qcdnum_params = QCDNUMParameters(order=2, α_S=0.118, q0=100.0, grid=grid,
                                 n_fixed_flav=5, iqc=1, iqb=1, iqt=1, weight_type=1);

# Initialise and set key parameters

QCDNUM.qcinit(-6, "")
QCDNUM.setord(qcdnum_params.order)
QCDNUM.setalf(qcdnum_params.α_S, qcdnum_params.q0)

# Build grids

g = qcdnum_params.grid
QCDNUM.gxmake(g.x_bounds, g.x_weights, g.x_num_bounds, g.nx,
              g.spline_interp);
QCDNUM.gqmake(g.qq_bounds, g.qq_weights, g.qq_num_bounds, g.nq);

# Define FFNS/VFNS

QCDNUM.setcbt(qcdnum_params.n_fixed_flav, qcdnum_params.iqc,
              qcdnum_params.iqb, qcdnum_params.iqt);

# Build weight tables

nw = QCDNUM.fillwt(qcdnum_params.weight_type)
nw = QCDNUM.zmfillw()

# ## Evolve the PDFs using QCDNUM
#
# Define input PDF function
# * See https://www.nikhef.nl/~h24/qcdnum-files/doc/qcdnum170115.pdf under `evolfg`
#
# There are functions available to generate the necessary input PDF function in
# the correct format for `QCDNUM.jl` (see `get_input_pdf_func()`), along with
# the mapping between this input function and quark species (see `input_pdf_map`).

# Get function and wrap with c-style pointer

my_func = get_input_pdf_func(hyper_params)
input_pdf = @cfunction(my_func, Float64, (Ref{Int32}, Ref{Float64}))

# Find index of starting scale and evolve

iq0 = QCDNUM.iqfrmq(qcdnum_params.q0)
pdf_loc = 1
eps = QCDNUM.evolfg(pdf_loc, input_pdf, input_pdf_map, iq0)

# ## Define necessary splines for cross section calculation

# For splines

splint_params = SPLINTParameters();
quark_coeffs = QuarkCoefficients();

# Define initial spline

QCDNUM.ssp_spinit(splint_params.nuser);
ia = QCDNUM.isp_s2make(splint_params.nsteps_x, splint_params.nsteps_q);
xnd = QCDNUM.ssp_unodes(ia, splint_params.nnodes_x, 0);
qnd = QCDNUM.ssp_vnodes(ia, splint_params.nnodes_q, 0);

# Check nodes and erase

QCDNUM.ssp_nprint(ia);
QCDNUM.ssp_erase(ia);

# Set nodes and fill spline with structure function

iaFLup = QCDNUM.isp_s2user(xnd, splint_params.nnodes_x, qnd, splint_params.nnodes_q);
QCDNUM.ssp_s2f123(iaFLup, pdf_loc, quark_coeffs.proup, 1, 0.0);

iaF2up = QCDNUM.isp_s2user(xnd, splint_params.nnodes_x, qnd, splint_params.nnodes_q);
QCDNUM.ssp_s2f123(iaF2up, pdf_loc, quark_coeffs.proup, 2, 0.0);

iaF3up = QCDNUM.isp_s2user(xnd, splint_params.nnodes_x, qnd, splint_params.nnodes_q);
QCDNUM.ssp_s2f123(iaF3up, pdf_loc, quark_coeffs.valup, 3, 0.0);
    
iaFLdn = QCDNUM.isp_s2user(xnd, splint_params.nnodes_x, qnd, splint_params.nnodes_q);
QCDNUM.ssp_s2f123(iaFLdn, pdf_loc, quark_coeffs.prodn, 1, 0.0);

iaF2dn = QCDNUM.isp_s2user(xnd, splint_params.nnodes_x, qnd, splint_params.nnodes_q);
QCDNUM.ssp_s2f123(iaF2dn, pdf_loc, quark_coeffs.prodn, 2, 0.0);

iaF3dn = QCDNUM.isp_s2user(xnd, splint_params.nnodes_x, qnd, splint_params.nnodes_q);
QCDNUM.ssp_s2f123(iaF3dn, 1, quark_coeffs.valdn, 3, 0.0);

# store spline addresses
QCDNUM.ssp_uwrite(splint_params.spline_addresses.F2up, Float64(iaF2up));
QCDNUM.ssp_uwrite(splint_params.spline_addresses.F2dn, Float64(iaF2dn));
QCDNUM.ssp_uwrite(splint_params.spline_addresses.F3up, Float64(iaF3up));
QCDNUM.ssp_uwrite(splint_params.spline_addresses.F3dn, Float64(iaF3dn));
QCDNUM.ssp_uwrite(splint_params.spline_addresses.FLup, Float64(iaFLup));
QCDNUM.ssp_uwrite(splint_params.spline_addresses.FLdn, Float64(iaFLdn));

my_func = get_input_xsec_func()
input_xsec = @cfunction(my_func, Float64, (Ref{Int32}, Ref{Int32}, Ref{UInt8}))

# plot

xsec_on_grid = zeros(g.nx, g.nq);

for ix = 1:g.nx
    for iq = 1:g.nq
        xsec_on_grid[ix, iq] = _fun_xsec_i(ix, iq)
    end
end

qcdnum_x_grid = QCDNUM.gxcopy(g.nx)
qcdnum_qq_grid =  QCDNUM.gqcopy(g.nq)
p1 = heatmap(qcdnum_x_grid, qcdnum_qq_grid, NaNMath.log10.(xsec_on_grid[:, :]'))
plot(p1, xlabel="x", ylabel="q2", 
    xaxis=:log, yaxis=:log)

plot(qcdnum_x_grid, NaNMath.log10.(xsec_on_grid[:, 4]), 
    label="Q2=141 (input scale)", lw=3)
plot!(qcdnum_x_grid, NaNMath.log10.(xsec_on_grid[:, 22]), label="Q2=1152", lw=3)
plot!(qcdnum_x_grid, NaNMath.log10.(xsec_on_grid[:, 35]), label="Q2=5233", lw=3)
plot!(qcdnum_x_grid, NaNMath.log10.(xsec_on_grid[:, 41]), label="Q2=10523", lw=3)
plot!(xaxis=:log, legend=:bottomleft, xlabel="x", 
    ylabel="log10(cross section spline input)", ylims=(-7, 5))

set_lepcharge(1)
iaF_eP = QCDNUM.isp_s2make(1, 2);
QCDNUM.ssp_uwrite(splint_params.spline_addresses.F_eP, Float64(iaF_eP));
QCDNUM.ssp_s2fill(iaF_eP, input_xsec, splint_params.rscut);

set_lepcharge(-1)
iaF_eM = QCDNUM.isp_s2make(1, 2);
QCDNUM.ssp_uwrite(splint_params.spline_addresses.F_eM, Float64(iaF_eM));
QCDNUM.ssp_s2fill(iaF_eM, input_xsec, splint_params.rscut);

# plot spline

spline = zeros(g.nx, g.nq);

for ix = 1:g.nx
    for iq = 1:g.nq
        spline[ix, iq] = QCDNUM.dsp_funs2(iaF_eP, qcdnum_x_grid[ix], 
            qcdnum_qq_grid[iq], 1)
    end
end

p1 = heatmap(qcdnum_x_grid, qcdnum_qq_grid, NaNMath.log10.(spline[:, :]'))
plot(p1, xlabel="x", ylabel="q2", 
    xaxis=:log, yaxis=:log)

## Integrate over the cross section spline and find expected events numbers

nbins = size(xbins_M_begin)[1]
IntXsec_eP = zeros(nbins);
IntXsec_eM = zeros(nbins);
for i in 1:nbins
    IntXsec_eP[i] = QCDNUM.dsp_ints2(iaF_eP, xbins_M_begin[i], xbins_M_end[i], 
        q2bins_M_begin[i], q2bins_M_end[i], 318., 4);
    IntXsec_eM[i] = QCDNUM.dsp_ints2(iaF_eM, xbins_M_begin[i], xbins_M_end[i], 
        q2bins_M_begin[i], q2bins_M_end[i], 318., 4);
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

xsec_pred_eM
