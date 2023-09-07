abstract type MetaData end

export MetaData


#function f_cross_section_to_counts(md::MetaData, integ_xsec_ep::Array{Float64},integ_xsec_em::Array{Float64},sys_err_params::Vector{Float64}=zeros(8)) 
#  println("D not call me!")
#end

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

function MetaDataZEUS(input::MetaDataZEUS, f::Float64) 
   output = input
   output.Ld_eMp= input.Ld_eMp*f
   output.Ld_ePp= input.Ld_ePp*f
   output
end 
