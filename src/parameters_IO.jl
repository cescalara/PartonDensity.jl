using HDF5
export pd_write_sim, pd_read_sim


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
