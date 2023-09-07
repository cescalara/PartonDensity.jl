export forward_model, forward_model_init

function forward_model(pdf_params::Union{BernsteinPDFParams,BernsteinDirichletPDFParams}, 
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
    nbins = size(m_xbins_M_begin)[1]
    integ_xsec_ep = zeros(nbins)
    integ_xsec_em = zeros(nbins)
    for i in 1:nbins
        integ_xsec_ep[i] = QCDNUM.dsp_ints2(iaF_eP, m_xbins_M_begin[i], m_xbins_M_end[i], m_q2bins_M_begin[i], m_q2bins_M_end[i], md.sqrtS, 4)
        integ_xsec_em[i] = QCDNUM.dsp_ints2(iaF_eM, m_xbins_M_begin[i], m_xbins_M_end[i], m_q2bins_M_begin[i], m_q2bins_M_end[i], md.sqrtS, 4)
    end

    counts_pred_ep, counts_pred_em = f_cross_section_to_counts(integ_xsec_ep,integ_xsec_em,sys_err_params)
    return counts_pred_ep, counts_pred_em
end
