export forward_model, forward_model_init

"""
    forward_model_init(qcdnum_grid, qcdnum_params)

Initialise forward model. Initialises QCDNUM and builds weight tables to 
save time in subsequent iterations. Must be called prior to first instance 
of `forward_model()`.
"""
function forward_model_init(qcdnum_params::QCDNUM.EvolutionParams, splint_params::QCDNUM.SPLINTParams)

    # Set up QCDNUM
    QCDNUM.init()
    QCDNUM.setord(qcdnum_params.order)
    QCDNUM.setalf(qcdnum_params.α_S, qcdnum_params.q0)

    # QCDNUM Grids
    nx, nq = QCDNUM.make_grid(qcdnum_params.grid_params)

    # Define FFNS/VFNS  
    QCDNUM.setcbt(qcdnum_params.n_fixed_flav, qcdnum_params.iqc, qcdnum_params.iqb, qcdnum_params.iqt)

    # Build weight tables
    nw = QCDNUM.fillwt(qcdnum_params.weight_type)
    nw = QCDNUM.zmfillw()

    # Define splines, nodes and addresses
    QCDNUM.ssp_spinit(splint_params.nuser)
    ia = QCDNUM.isp_s2make(splint_params.nsteps_x, splint_params.nsteps_q)
    xnd = QCDNUM.ssp_unodes(ia, splint_params.nnodes_x, 0)
    qnd = QCDNUM.ssp_vnodes(ia, splint_params.nnodes_q, 0)
    QCDNUM.ssp_erase(ia)

    iaFLup = QCDNUM.isp_s2user(xnd, splint_params.nnodes_x, qnd, splint_params.nnodes_q)
    iaF2up = QCDNUM.isp_s2user(xnd, splint_params.nnodes_x, qnd, splint_params.nnodes_q)
    iaF3up = QCDNUM.isp_s2user(xnd, splint_params.nnodes_x, qnd, splint_params.nnodes_q)
    iaFLdn = QCDNUM.isp_s2user(xnd, splint_params.nnodes_x, qnd, splint_params.nnodes_q)
    iaF2dn = QCDNUM.isp_s2user(xnd, splint_params.nnodes_x, qnd, splint_params.nnodes_q)
    iaF3dn = QCDNUM.isp_s2user(xnd, splint_params.nnodes_x, qnd, splint_params.nnodes_q)
    iaF_eP = QCDNUM.isp_s2make(1, 2)
    iaF_eM = QCDNUM.isp_s2make(1, 2)

    # Store spline addresses
    QCDNUM.ssp_uwrite(splint_params.spline_addresses.F2up, Float64(iaF2up))
    QCDNUM.ssp_uwrite(splint_params.spline_addresses.F2dn, Float64(iaF2dn))
    QCDNUM.ssp_uwrite(splint_params.spline_addresses.F3up, Float64(iaF3up))
    QCDNUM.ssp_uwrite(splint_params.spline_addresses.F3dn, Float64(iaF3dn))
    QCDNUM.ssp_uwrite(splint_params.spline_addresses.FLup, Float64(iaFLup))
    QCDNUM.ssp_uwrite(splint_params.spline_addresses.FLdn, Float64(iaFLdn))
    QCDNUM.ssp_uwrite(splint_params.spline_addresses.F_eP, Float64(iaF_eP))
    QCDNUM.ssp_uwrite(splint_params.spline_addresses.F_eM, Float64(iaF_eM))

end


"""
    forward_model(pdf_params, qcdnum_grid, splint_params, quark_coeffs, sys_err_params)

Go from input PDF parameters to the expected number of events in bins.
"""
function forward_model(pdf_params::AbstractPDFParams,
                       qcdnum_params::QCDNUM.EvolutionParams,
                       splint_params::QCDNUM.SPLINTParams, 
                       quark_coeffs::QuarkCoefficients,
                       md::MetaData = MD_ZEUS_I1787035
                       ,sys_err_params::Vector{Float64}=zeros(nsyst))

    # Get input PDF function
    my_func = get_input_pdf_func(pdf_params)
    input_pdf = QCDNUM.InputPDF(func=my_func, map=input_pdf_map)

    # Evolve PDFs
    iq0 = QCDNUM.iqfrmq(qcdnum_params.q0)
    ϵ = QCDNUM.evolfg(qcdnum_params.output_pdf_loc, input_pdf.cfunc, input_pdf.map, iq0)

    # Debugging
    if ϵ > 0.05

        @warn "QCDNUM.evolfg(): Spline issues detected" ϵ pdf_params

    end

    # Read spline addresses
    iaF2up = Int64(QCDNUM.dsp_uread(splint_params.spline_addresses.F2up))
    iaF2dn = Int64(QCDNUM.dsp_uread(splint_params.spline_addresses.F2dn))
    iaF3up = Int64(QCDNUM.dsp_uread(splint_params.spline_addresses.F3up))
    iaF3dn = Int64(QCDNUM.dsp_uread(splint_params.spline_addresses.F3dn))
    iaFLup = Int64(QCDNUM.dsp_uread(splint_params.spline_addresses.FLup))
    iaFLdn = Int64(QCDNUM.dsp_uread(splint_params.spline_addresses.FLdn))
    iaF_eP = Int64(QCDNUM.dsp_uread(splint_params.spline_addresses.F_eP))
    iaF_eM = Int64(QCDNUM.dsp_uread(splint_params.spline_addresses.F_eM))

    # Splines for structure functions
    QCDNUM.ssp_s2f123(iaFLup, qcdnum_params.output_pdf_loc, quark_coeffs.proup, 1, 0.0)
    QCDNUM.ssp_s2f123(iaF2up, qcdnum_params.output_pdf_loc, quark_coeffs.proup, 2, 0.0)
    QCDNUM.ssp_s2f123(iaF3up, qcdnum_params.output_pdf_loc, quark_coeffs.valup, 3, 0.0)
    QCDNUM.ssp_s2f123(iaFLdn, qcdnum_params.output_pdf_loc, quark_coeffs.prodn, 1, 0.0)
    QCDNUM.ssp_s2f123(iaF2dn, qcdnum_params.output_pdf_loc, quark_coeffs.prodn, 2, 0.0)
    QCDNUM.ssp_s2f123(iaF3dn, qcdnum_params.output_pdf_loc, quark_coeffs.valdn, 3, 0.0)

    # Get input cross section function
    my_funcp = get_input_xsec_func(1, md)
    input_xsecp = @cfunction($my_funcp, Float64, (Ref{Int32}, Ref{Int32}, Ref{UInt8}))

    my_funcm = get_input_xsec_func(-1, md)
    input_xsecm = @cfunction($my_funcm, Float64, (Ref{Int32}, Ref{Int32}, Ref{UInt8}))

    # Make two cross section splines
    QCDNUM.ssp_s2fill(iaF_eP, input_xsecp, splint_params.rscut)

    QCDNUM.ssp_s2fill(iaF_eM, input_xsecm, splint_params.rscut)

    # Integrate over cross section
    nbins = size(md.m_xbins_M_begin)[1]
    bins_axis = axes(md.m_xbins_M_begin, 1)

    integ_xsec_ep = similar(md.m_xbins_M_begin)
    integ_xsec_em = similar(md.m_xbins_M_begin)
    for i in bins_axis
        integ_xsec_ep[i] = QCDNUM.dsp_ints2(iaF_eP, md.m_xbins_M_begin[i], md.m_xbins_M_end[i], md.m_q2bins_M_begin[i], md.m_q2bins_M_end[i], md.sqrtS, 4)
        integ_xsec_em[i] = QCDNUM.dsp_ints2(iaF_eM, md.m_xbins_M_begin[i], md.m_xbins_M_end[i], md.m_q2bins_M_begin[i], md.m_q2bins_M_end[i], md.sqrtS, 4)
    end

    counts_pred_ep, counts_pred_em = f_cross_section_to_counts(integ_xsec_ep,integ_xsec_em,sys_err_params)
    return counts_pred_ep, counts_pred_em
end
