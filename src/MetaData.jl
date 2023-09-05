export MetaData
export PDData
export MD_DUMMY
export PD_DUMMY

struct PDData
q2bins_M_begin::Array{Float64}
q2bins_M_end::Array{Float64}
xbins_M_begin::Array{Float64}
xbins_M_end::Array{Float64}
end

const PD_DUMMY = PDData(Array{Float64,1}(undef, 1),Array{Float64,1}(undef, 1),Array{Float64,1}(undef, 1),Array{Float64,1}(undef, 1));

struct MetaData
 name::String
 Ld_ePp :: Float64
 Ld_eMp :: Float64
 Ld_ePp_uncertainty :: Float64
 Ld_eMp_uncertainty :: Float64
 sqrtS::Float64
 nsyst::Int
 D::PDData
end

const MD_DUMMY = MetaData("DUMMY", 141.44, 185.018, 0.018,0.018, 318.0, 8,PD_DUMMY);

