# # Forward model
# 
# Here, we go through an example of simulating the full forward model,
# from the prior definition to the expected number of events in different
# bins of the detector response.

using QCDNUM, PartonDensity
using Plots, Printf, NaNMath, Parameters, Random, Distributions

zeus_include_path = string(chop(pathof(PartonDensity), tail=20), "data/ZEUS_I1787035/ZEUS_I1787035.jl")
MD_ZEUS_I1787035 = include(zeus_include_path)

const MD_DOCS = MD_ZEUS_I1787035

# ## Define input PDFs
# We can use `DirichletPDFParams` or `ValencePDFParams`, as long
# as we do so according to the *PDF parametrisation and priors* docs.
random_seed = 42

weights = [3.0, 1.0, 5.0, 5.0, 1.0, 1.0, 1.0, 0.5, 0.5]
θ = rand(MersenneTwister(random_seed), Dirichlet(weights))
pdf_params = DirichletPDFParams(K_u=4.0, K_d=6.0, λ_g1=0.7, λ_g2=-0.4,
    K_g=6.0, λ_q=-0.5, K_q=5.0, θ=θ);

# Plot the input PDFs

plot_input_pdfs(pdf_params)

# Sanity check that sum = 1

int_xtotx(pdf_params) ≈ 1

# ## Define QCDNUM grids, weights and settings

grid = QCDNUM.GridParams(x_min=[1.0e-3], x_weights=[1], nx=100,
    qq_bounds=[1.0e2, 3.0e4], qq_weights=[1.0, 1.0], nq=50, spline_interp=3)

qcdnum_params = QCDNUM.EvolutionParams(order=2, α_S=0.118, q0=100.0, grid_params=grid,
    n_fixed_flav=5, iqc=1, iqb=1, iqt=1, weight_type=1);

# Initialise

QCDNUM.init()

# ## Evolve the PDFs using QCDNUM
#
# Define input PDF function
# * See QCDNUM docs under `evolfg`
#
# There are functions available to generate the necessary input PDF function in
# the correct format for `QCDNUM.jl` (see `get_input_pdf_func()`), along with
# the mapping between this input function and quark species (see `input_pdf_map`).

# Get function and PDF input map to fully describe the `QCDNUM.InputPDF`

my_func = get_input_pdf_func(pdf_params)
input_pdf = QCDNUM.InputPDF(func=my_func, map=input_pdf_map)

# Evolve the PDF over the specified grid

ϵ = QCDNUM.evolve(input_pdf, qcdnum_params)

# For splines
nw = QCDNUM.zmfillw()
splint_params = QCDNUM.SPLINTParams();
quark_coeffs = QuarkCoefficients();

# Define initial spline

QCDNUM.ssp_spinit(splint_params.nuser)

ia = QCDNUM.isp_s2make(splint_params.nsteps_x, splint_params.nsteps_q);
xnd = QCDNUM.ssp_unodes(ia, splint_params.nnodes_x, 0);
qnd = QCDNUM.ssp_vnodes(ia, splint_params.nnodes_q, 0);

# Check nodes and erase

QCDNUM.ssp_nprint(ia);
QCDNUM.ssp_erase(ia);

# Set nodes and fill spline with structure function

iaFLup = QCDNUM.isp_s2user(xnd, splint_params.nnodes_x, qnd, splint_params.nnodes_q);
QCDNUM.ssp_s2f123(iaFLup, qcdnum_params.output_pdf_loc, quark_coeffs.proup, 1, 0.0);

iaF2up = QCDNUM.isp_s2user(xnd, splint_params.nnodes_x, qnd, splint_params.nnodes_q);
QCDNUM.ssp_s2f123(iaF2up, qcdnum_params.output_pdf_loc, quark_coeffs.proup, 2, 0.0);

iaF3up = QCDNUM.isp_s2user(xnd, splint_params.nnodes_x, qnd, splint_params.nnodes_q);
QCDNUM.ssp_s2f123(iaF3up, qcdnum_params.output_pdf_loc, quark_coeffs.valup, 3, 0.0);

iaFLdn = QCDNUM.isp_s2user(xnd, splint_params.nnodes_x, qnd, splint_params.nnodes_q);
QCDNUM.ssp_s2f123(iaFLdn, qcdnum_params.output_pdf_loc, quark_coeffs.prodn, 1, 0.0);

iaF2dn = QCDNUM.isp_s2user(xnd, splint_params.nnodes_x, qnd, splint_params.nnodes_q);
QCDNUM.ssp_s2f123(iaF2dn, qcdnum_params.output_pdf_loc, quark_coeffs.prodn, 2, 0.0);

iaF3dn = QCDNUM.isp_s2user(xnd, splint_params.nnodes_x, qnd, splint_params.nnodes_q);
QCDNUM.ssp_s2f123(iaF3dn, qcdnum_params.output_pdf_loc, quark_coeffs.valdn, 3, 0.0);

# store spline addresses
QCDNUM.ssp_uwrite(splint_params.spline_addresses.F2up, Float64(iaF2up));
QCDNUM.ssp_uwrite(splint_params.spline_addresses.F2dn, Float64(iaF2dn));
QCDNUM.ssp_uwrite(splint_params.spline_addresses.F3up, Float64(iaF3up));
QCDNUM.ssp_uwrite(splint_params.spline_addresses.F3dn, Float64(iaF3dn));
QCDNUM.ssp_uwrite(splint_params.spline_addresses.FLup, Float64(iaFLup));
QCDNUM.ssp_uwrite(splint_params.spline_addresses.FLdn, Float64(iaFLdn));

my_funcp = get_input_xsec_func(1, MD_DOCS) # charge = 1
input_xsecp = @cfunction(my_funcp, Float64, (Ref{Int32}, Ref{Int32}, Ref{UInt8}))

my_funcm = get_input_xsec_func(-1, MD_DOCS) # charge = -1
input_xsecm = @cfunction(my_funcm, Float64, (Ref{Int32}, Ref{Int32}, Ref{UInt8}))

# plot
g = qcdnum_params.grid_params
xsec_on_grid = zeros(g.nx, g.nq);

for ix = 1:g.nx
    for iq = 1:g.nq
        xsec_on_grid[ix, iq] = _fun_xsec_i(1, MD_DOCS, ix, iq) # charge = 1
    end
end

qcdnum_x_grid = QCDNUM.gxcopy(g.nx)
qcdnum_qq_grid = QCDNUM.gqcopy(g.nq)
p1 = heatmap(qcdnum_x_grid, qcdnum_qq_grid, NaNMath.log10.(xsec_on_grid[:, :]'))
plot(p1, xlabel="x", ylabel="q2",
    xaxis=:log, yaxis=:log)

#

plot(qcdnum_x_grid, NaNMath.log10.(xsec_on_grid[:, 4]),
    label="Q2=141 (input scale)", lw=3)
plot!(qcdnum_x_grid, NaNMath.log10.(xsec_on_grid[:, 22]), label="Q2=1152", lw=3)
plot!(qcdnum_x_grid, NaNMath.log10.(xsec_on_grid[:, 35]), label="Q2=5233", lw=3)
plot!(qcdnum_x_grid, NaNMath.log10.(xsec_on_grid[:, 41]), label="Q2=10523", lw=3)
plot!(xaxis=:log, legend=:bottomleft, xlabel="x",
    ylabel="log10(cross section spline input)", ylims=(-7, 5))

#

iaF_eP = QCDNUM.isp_s2make(1, 2);
QCDNUM.ssp_uwrite(splint_params.spline_addresses.F_eP, Float64(iaF_eP));
QCDNUM.ssp_s2fill(iaF_eP, input_xsecp, splint_params.rscut);

iaF_eM = QCDNUM.isp_s2make(1, 2);
QCDNUM.ssp_uwrite(splint_params.spline_addresses.F_eM, Float64(iaF_eM));
QCDNUM.ssp_s2fill(iaF_eM, input_xsecm, splint_params.rscut);

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

# ## Integrate over the cross section spline and find expected events numbers
#
# Here, we neglect any possible contribution from systematic errors

nbins = size(MD_DOCS.m_xbins_M_begin)[1]
IntXsec_eP = zeros(nbins);
IntXsec_eM = zeros(nbins);
for i in 1:nbins
    IntXsec_eP[i] = QCDNUM.dsp_ints2(iaF_eP, MD_DOCS.m_xbins_M_begin[i], MD_DOCS.m_xbins_M_end[i], MD_DOCS.m_q2bins_M_begin[i], MD_DOCS.m_q2bins_M_end[i], MD_DOCS.sqrtS, 4)
    IntXsec_eM[i] = QCDNUM.dsp_ints2(iaF_eM, MD_DOCS.m_xbins_M_begin[i], MD_DOCS.m_xbins_M_end[i], MD_DOCS.m_q2bins_M_begin[i], MD_DOCS.m_q2bins_M_end[i], MD_DOCS.sqrtS, 4)
end

counts_pred_eP, counts_pred_eM = MD_DOCS.f_cross_section_to_counts(MD_DOCS.Ld_ePp, MD_DOCS.Ld_eMp, IntXsec_eP, IntXsec_eM)

counts_pred_eM
