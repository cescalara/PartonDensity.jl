export forward_model, forward_model_init

"""
    forward_model_init(qcdnum_grid, qcdnum_params)

Initialise forward model. Initialises QCDNUM and builds weight tables to 
save time in subsequent iterations. Must be called prior to first instance 
of `forward_model()`.
"""
function forward_model_init(qcdnum_grid::QCDNUMGrid, qcdnum_params::QCDNUMParameters)

    # Set up QCDNUM
    QCDNUM.qcinit(-6, "")
    QCDNUM.setord(qcdnum_params.order)
    QCDNUM.setalf(qcdnum_params.α_S, qcdnum_params.q0)

    # QCDNUM Grids
    g = qcdnum_params.grid
    QCDNUM.gxmake(g.x_min, g.x_weights, g.x_num_bounds, g.nx,
                  g.spline_interp)
    QCDNUM.gqmake(g.qq_bounds, g.qq_weights, g.qq_num_bounds, g.nq)

    # Define FFNS/VFNS
    QCDNUM.setcbt(qcdnum_params.n_fixed_flav, qcdnum_params.iqc,
                  qcdnum_params.iqb, qcdnum_params.iqt)

    # Build weight tables
    # TODO: Use saved weights file once QCXDNUM fixed
    nw = QCDNUM.fillwt(qcdnum_params.weight_type)
    nw = QCDNUM.zmfillw()

end


"""
    forward_model(pdf_params, qcdnum_grid, 
                  splint_params, quark_coeffs)

Go from input PDF parameters to the expected number of events in bins.
"""
function forward_model(pdf_params::PDFParameters, qcdnum_params::QCDNUMParameters,
                       splint_params::SPLINTParameters, quark_coeffs::QuarkCoefficients)


    # Get input PDF function
    my_func = get_input_pdf_func(pdf_params)
    input_pdf = @cfunction($my_func, Float64, (Ref{Int32}, Ref{Float64}))

    # Evolve PDFs
    iq0 = QCDNUM.iqfrmq(qcdnum_params.q0)
    pdf_loc = 1
    eps = QCDNUM.evolfg(pdf_loc, input_pdf, input_pdf_map, iq0)

    # Define splines for structure functions
    QCDNUM.ssp_spinit(splint_params.nuser)
    ia = QCDNUM.isp_s2make(splint_params.nsteps_x, splint_params.nsteps_q)
    xnd = QCDNUM.ssp_unodes(ia, splint_params.nnodes_x, 0)
    qnd = QCDNUM.ssp_vnodes(ia, splint_params.nnodes_q, 0)
    QCDNUM.ssp_erase(ia);

    iaFLup = QCDNUM.isp_s2user(xnd, splint_params.nnodes_x, qnd, splint_params.nnodes_q)
    QCDNUM.ssp_s2f123(iaFLup, pdf_loc, quark_coeffs.proup, 1, 0.0)

    iaF2up = QCDNUM.isp_s2user(xnd, splint_params.nnodes_x, qnd, splint_params.nnodes_q)
    QCDNUM.ssp_s2f123(iaF2up, pdf_loc, quark_coeffs.proup, 2, 0.0)

    iaF3up = QCDNUM.isp_s2user(xnd, splint_params.nnodes_x, qnd, splint_params.nnodes_q)
    QCDNUM.ssp_s2f123(iaF3up, pdf_loc, quark_coeffs.valup, 3, 0.0)
    
    iaFLdn = QCDNUM.isp_s2user(xnd, splint_params.nnodes_x, qnd, splint_params.nnodes_q)
    QCDNUM.ssp_s2f123(iaFLdn, pdf_loc, quark_coeffs.prodn, 1, 0.0)
    
    iaF2dn = QCDNUM.isp_s2user(xnd, splint_params.nnodes_x, qnd, splint_params.nnodes_q)
    QCDNUM.ssp_s2f123(iaF2dn, pdf_loc, quark_coeffs.prodn, 2, 0.0)

    iaF3dn = QCDNUM.isp_s2user(xnd, splint_params.nnodes_x, qnd, splint_params.nnodes_q)
    QCDNUM.ssp_s2f123(iaF3dn, 1, quark_coeffs.valdn, 3, 0.0)

    # Store spline addresses
    QCDNUM.ssp_uwrite(splint_params.spline_addresses.F2up, Float64(iaF2up))
    QCDNUM.ssp_uwrite(splint_params.spline_addresses.F2dn, Float64(iaF2dn))
    QCDNUM.ssp_uwrite(splint_params.spline_addresses.F3up, Float64(iaF3up))
    QCDNUM.ssp_uwrite(splint_params.spline_addresses.F3dn, Float64(iaF3dn))
    QCDNUM.ssp_uwrite(splint_params.spline_addresses.FLup, Float64(iaFLup))
    QCDNUM.ssp_uwrite(splint_params.spline_addresses.FLdn, Float64(iaFLdn))

    # Get input cross section function
    my_func = get_input_xsec_func()
    input_xsec = @cfunction($my_func, Float64, (Ref{Int32}, Ref{Int32}, Ref{UInt8}))

    # Make two cross section splines
    set_lepcharge(1)
    iaF_eP = QCDNUM.isp_s2make(1, 2)
    QCDNUM.ssp_uwrite(splint_params.spline_addresses.F_eP, Float64(iaF_eP))
    QCDNUM.ssp_s2fill(iaF_eP, input_xsec, splint_params.rscut)

    set_lepcharge(-1)
    iaF_eM = QCDNUM.isp_s2make(1, 2)
    QCDNUM.ssp_uwrite(splint_params.spline_addresses.F_eM, Float64(iaF_eM))
    QCDNUM.ssp_s2fill(iaF_eM, input_xsec, splint_params.rscut)

    # Integrate over cross section
    nbins = size(xbins_M_begin)[1]
    integ_xsec_ep = zeros(nbins);
    integ_xsec_em = zeros(nbins);

    for i in 1:nbins
        integ_xsec_ep[i] = QCDNUM.dsp_ints2(iaF_eP, xbins_M_begin[i], xbins_M_end[i],
                                            q2bins_M_begin[i], q2bins_M_end[i], 318., 4)
        integ_xsec_em[i] = QCDNUM.dsp_ints2(iaF_eM, xbins_M_begin[i], xbins_M_end[i],
                                            q2bins_M_begin[i], q2bins_M_end[i], 318., 4)
    end

    # Fold through response to get counts
    ePp = 0
    eMp = 1

    TM_eP = get_TM_elements(ePp)
    TM_eM = get_TM_elements(eMp)

    K_eP = get_K_elements(ePp)
    K_eM = get_K_elements(eMp)

    nbins_out = size(TM_eP)[2]

    counts_pred_ep = zeros(nbins_out)
    counts_pred_em = zeros(nbins_out)

    for j in 1:nbins_out
        
        for i in 1:nbins

            counts_pred_ep[j] += TM_eP[i, j] * (1.0/K_eP[i]) * integ_xsec_ep[i]
            counts_pred_em[j] += TM_eM[i, j] * (1.0/K_eM[i]) * integ_xsec_em[i]
            
        end

    end

    return counts_pred_ep, counts_pred_em
    
end
