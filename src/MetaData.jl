abstract type MetaData end

export MetaData
export scale_lumi

mutable struct MetaDataZEUS{T} <: MetaData
 name::String
 Ld_ePp :: Float64
 Ld_eMp :: Float64
 Ld_ePp_uncertainty :: Float64
 Ld_eMp_uncertainty :: Float64
 sqrtS::Float64
 m_q2bins_M_begin::Array{Float64}
 m_q2bins_M_end::Array{Float64}
 m_xbins_M_begin::Array{Float64}
 m_xbins_M_end::Array{Float64}
 m_Data_Events_eMp::Array{Int64}
 m_Data_Events_ePp::Array{Int64}
 f_cross_section_to_counts::T
end

mutable struct MetaDataIO <: MetaData
 name::String
 Ld_ePp :: Float64
 Ld_eMp :: Float64
 Ld_ePp_uncertainty :: Float64
 Ld_eMp_uncertainty :: Float64
 sqrtS::Float64
end

function MetaDataIO(input::MetaData) 
   output::MetaDataIO  = MetaDataIO(input.name,
   input.Ld_eMp,
   input.Ld_ePp,
   input.Ld_eMp_uncertainty,
   input.Ld_ePp_uncertainty,
   input.sqrtS
)   
   output
end 


function scale_lumi(input, f::Float64) 
   output = input
   output.Ld_eMp= input.Ld_eMp*f
   output.Ld_ePp= input.Ld_ePp*f
   output
end 
