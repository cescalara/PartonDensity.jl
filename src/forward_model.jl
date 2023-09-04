using HDF5

export forward_model, forward_model_init
export pd_write_sim, pd_read_sim


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
    QCDNUM.setcbt(qcdnum_params.n_fixed_flav, qcdnum_params.iqc,
        qcdnum_params.iqb, qcdnum_params.iqt)

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
    forward_model(pdf_params, qcdnum_grid, 
                  splint_params, quark_coeffs, sys_err_params)

Go from input PDF parameters to the expected number of events in bins.
"""
function forward_model(pdf_params::AbstractPDFParams, 
                       qcdnum_params::QCDNUM.EvolutionParams,
                       splint_params::QCDNUM.SPLINTParams, 
                       quark_coeffs::QuarkCoefficients,
                       md::MetaData = MD_ZEUS,
                       sys_err_params::Vector{Float64}=zeros(nsyst))

    # Get input PDF function
    my_func = get_input_pdf_func(pdf_params)
    input_pdf = QCDNUM.InputPDF(func=my_func, map=input_pdf_map)

    # Evolve PDFs
    iq0 = QCDNUM.iqfrmq(qcdnum_params.q0)

    ϵ = QCDNUM.evolfg(qcdnum_params.output_pdf_loc, input_pdf.cfunc, input_pdf.map, iq0)

    # Debugging
    if ϵ > 0.05

        @warn "QCDNUM.evolfg(): Spline issues detected" eps pdf_params

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
    QCDNUM.ssp_s2f123(iaF3dn, qcdnum_params.output_pdf a dozen of trivial MRs_loc, quark_coeffs.valdn, 3, 0.0)

    # Get input cross section function
    my_funcp = get_input_xsec_func(1)
    input_xsecp = @cfunction($my_funcp, Float64, (Ref{Int32}, Ref{Int32}, Ref{UInt8}))

    my_funcm = get_input_xsec_func(-1)
    input_xsecm = @cfunction($my_funcm, Float64, (Ref{Int32}, Ref{Int32}, Ref{UInt8}))

    # Make two cross section splines
    QCDNUM.ssp_s2fill(iaF_eP, input_xsecp, splint_params.rscut)

    QCDNUM.ssp_s2fill(iaF_eM, input_xsecm, splint_params.rscut)

    # Integrate over cross section
    nbins = size(xbins_M_begin)[1]
    bins_axis = axes(xbins_M_begin, 1)

    integ_xsec_ep = similar(xbins_M_begin)
    integ_xsec_em = similar(xbins_M_begin)
    for i in bins_axis
        integ_xsec_ep[i] = QCDNUM.dsp_ints2(iaF_eP, xbins_M_begin[i], xbins_M_end[i], q2bins_M_begin[i], q2bins_M_end[i], md.sqrtS, 4)
        integ_xsec_em[i] = QCDNUM.dsp_ints2(iaF_eM, xbins_M_begin[i], xbins_M_end[i], q2bins_M_begin[i], q2bins_M_end[i], md.sqrtS, 4)
    end

    # Fold through response to get counts
    ePp = 0
    eMp = 1

    TM_eP = get_TM_elements(ePp,md)
    TM_eM = get_TM_elements(eMp,md)


    K_eP = get_K_elements(ePp)
    K_eM = get_K_elements(eMp)

    T = promote_type(map(eltype, (sys_err_params, Tnm_sys_ePp, Tnm_sys_eMp,TM_eP, TM_eM, K_eP, K_eM, integ_xsec_ep, integ_xsec_em))...)

    counts_pred_ep = similar(TM_eP, T, size(TM_eP, 2))
    counts_pred_em = similar(TM_eM, T, size(TM_eM, 2))

    syserr_axis = axes(sys_err_params, 1)

    @argcheck axes(TM_eP, 2) == axes(TM_eM, 2) == axes(Tnm_sys_ePp, 2) == axes(Tnm_sys_eMp, 2)
    @argcheck axes(TM_eP, 1) == axes(TM_eM, 1) == axes(K_eP, 1) == axes(K_eM, 1)
    @argcheck axes(sys_err_params, 1) == axes(Tnm_sys_ePp, 3) == axes(Tnm_sys_eMp, 3)

    bin_out_axis = axes(counts_pred_ep, 1)
    bin_axis = axes(TM_eP, 1)

    # Calculate TotSys_var_em == sys_err_params[k] * Tnm_sys_ePp[i,j,k] up front?

    fill!(counts_pred_ep, 0)
    fill!(counts_pred_em, 0)

    for j in bin_out_axis
        TotSys_var_ep::T = 0
        TotSys_var_em::T = 0
        for i in bin_axis
            # Add variation for parameters, will do nothing if no systematic errors are provided:
            for k in syserr_axis
                TotSys_var_ep += sys_err_params[k] * Tnm_sys_ePp[i, j, k]
                TotSys_var_em += sys_err_params[k] * Tnm_sys_eMp[i, j, k]
            end
            counts_pred_ep[j] += (TM_eP[i, j] + TotSys_var_ep) * (1.0 / K_eP[i]) * integ_xsec_ep[i]
            counts_pred_em[j] += (TM_eM[i, j] + TotSys_var_em) * (1.0 / K_eM[i]) * integ_xsec_em[i]
        end
    end

    return counts_pred_ep, counts_pred_em
end

"""
    pd_write_sim(file_name, pdf_params, sim_data)

Store the simulation truth and simulated data in an HDF5 file.
"""
function pd_write_sim(file_name::String, pdf_params::Union{ValencePDFParams,DirichletPDFParams}, sim_data::Dict{String,Any})

    h5open(file_name, "w") do fid

        # store sim_data
        sim_data_group = create_group(fid, "data")
        for (key, value) in sim_data
            sim_data_group[key] = value
        end

        # store pdf_params
        truth_group = create_group(fid, "truth")
        truth_group["λ_u"] = pdf_params.λ_u
        truth_group["K_u"] = pdf_params.K_u
        truth_group["λ_d"] = pdf_params.λ_d
        truth_group["K_d"] = pdf_params.K_d
        truth_group["λ_g1"] = pdf_params.λ_g1
        truth_group["λ_g2"] = pdf_params.λ_g2
        truth_group["K_g"] = pdf_params.K_g
        truth_group["λ_q"] = pdf_params.λ_q
        truth_group["K_q"] = pdf_params.K_q
        truth_group["θ"] = pdf_params.θ
        truth_group["param_type"] = pdf_params.param_type

    end

    return nothing
end

"""
    pd_read_sim(file_name)

Read in the simulated truth and simulated data from HDF5 file.
"""
function pd_read_sim(file_name::String)

    local pdf_params
    sim_data = Dict{String,Any}()

    h5open(file_name, "r") do fid

        # read sim_data
        for (key, value) in zip(keys(fid["data"]), fid["data"])
            sim_data[key] = read(value)
        end

        g = fid["truth"]
        if read(g["param_type"]) == VALENCE_TYPE

            pdf_params = ValencePDFParams(λ_u=read(g["λ_u"]), K_u=read(g["K_u"]),
                λ_d=read(g["λ_d"]), K_d=read(g["K_d"]),
                λ_g1=read(g["λ_g1"]), λ_g2=read(g["λ_g2"]),
                K_g=read(g["K_g"]), λ_q=read(g["λ_q"]), K_q=read(g["K_q"]),
                θ=read(g["θ"]))

        elseif read(g["param_type"]) == DIRICHLET_TYPE

            pdf_params = DirichletPDFParams(λ_u=read(g["λ_u"]), K_u=read(g["K_u"]),
                λ_d=read(g["λ_d"]), K_d=read(g["K_d"]),
                λ_g1=read(g["λ_g1"]), λ_g2=read(g["λ_g2"]),
                K_g=read(g["K_g"]), λ_q=read(g["λ_q"]), K_q=read(g["K_q"]),
                θ=read(g["θ"]))

        elseif read(g["param_type"]) == BERNSTEIN_TYPE

            pdf_params = BernsteinPDFParams(U_list=read(g["U_list"]),
                D_list=read(g["D_list"]),
                λ_g1=read(g["λ_g1"]), λ_g2=read(g["λ_g2"]),
                K_g=read(g["K_g"]), λ_q=read(g["λ_q"]), K_q=read(g["K_q"]),
                weights=read(g["weights"]),
                θ=read(g["θ"]))

        elseif read(g["param_type"]) == BERNSTEIN_DIRICHLET_TYPE

            pdf_params = BernsteinDirichletPDFParams(U_list=read(g["U_list"]),
                D_list=read(g["D_list"]),
                λ_g1=read(g["λ_g1"]), λ_g2=read(g["λ_g2"]),
                K_g=read(g["K_g"]), λ_q=read(g["λ_q"]), K_q=read(g["K_q"]),
                weights=read(g["weights"]),
                θ=read(g["θ"]))

        else

            @error "PDF parametrisation not recognised."

        end

    end

    return pdf_params, sim_data

end
