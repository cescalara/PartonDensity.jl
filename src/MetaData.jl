export MetaData

struct MetaData
 name::String
 Ld_ePp :: Float64
 Ld_eMp :: Float64
 Ld_ePp_uncertainty :: Float64
 Ld_eMp_uncertainty :: Float64
 sqrtS::Float64
end

function MetaData(input::MetaData, f::Float64) 
   output = input
   output.Ld_eMp= input.Ld_eMp*f
   output.Ld_ePp= input.Ld_ePp*f
   output
end 