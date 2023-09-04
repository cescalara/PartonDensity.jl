export SysPenality, Init_sys
export Tnm_sys_ePp, Tnm_sys_eMp
export MD_ZEUS

include("ePp_jl/SysTnm_Eehigh_ePp.jl")
include("ePp_jl/SysTnm_Eelow_ePp.jl")
include("ePp_jl/SysTnm_Eeconehigh_ePp.jl")
include("ePp_jl/SysTnm_Eeconelow_ePp.jl")
include("ePp_jl/SysTnm_Eereshigh_ePp.jl")
include("ePp_jl/SysTnm_Eereslow_ePp.jl")
include("ePp_jl/SysTnm_Ejhigh_ePp.jl")
include("ePp_jl/SysTnm_Ejlow_ePp.jl")
include("ePp_jl/SysTnm_FBcalhigh_ePp.jl")
include("ePp_jl/SysTnm_FBcallow_ePp.jl")
include("ePp_jl/SysTnm_FCALxhigh_ePp.jl")
include("ePp_jl/SysTnm_FCALxlow_ePp.jl")
include("ePp_jl/SysTnm_FCALyhigh_ePp.jl")
include("ePp_jl/SysTnm_FCALylow_ePp.jl")
include("ePp_jl/SysTnm_AriMepsUp_ePp.jl")
include("ePp_jl/SysTnm_AriMepsDown_ePp.jl")

include("eMp_jl/SysTnm_Eehigh_eMp.jl")
include("eMp_jl/SysTnm_Eelow_eMp.jl")
include("eMp_jl/SysTnm_Eeconehigh_eMp.jl")
include("eMp_jl/SysTnm_Eeconelow_eMp.jl")
include("eMp_jl/SysTnm_Eereshigh_eMp.jl")
include("eMp_jl/SysTnm_Eereslow_eMp.jl")
include("eMp_jl/SysTnm_Ejhigh_eMp.jl")
include("eMp_jl/SysTnm_Ejlow_eMp.jl")
include("eMp_jl/SysTnm_FBcalhigh_eMp.jl")
include("eMp_jl/SysTnm_FBcallow_eMp.jl")
include("eMp_jl/SysTnm_FCALxhigh_eMp.jl")
include("eMp_jl/SysTnm_FCALxlow_eMp.jl")
include("eMp_jl/SysTnm_FCALyhigh_eMp.jl")
include("eMp_jl/SysTnm_FCALylow_eMp.jl")
include("eMp_jl/SysTnm_AriMepsUp_eMp.jl")
include("eMp_jl/SysTnm_AriMepsDown_eMp.jl")
# include("../../zeus_bin_edges.jl")


const NP1 = size(Tnm_Eehigh_ePp, 1)
const NP2 = size(Tnm_Eehigh_ePp, 2)
const nsyst = 8


const Tnm_sys_ePp = zeros(NP1, NP2, nsyst)
const Tnm_sys_eMp = zeros(NP1, NP2, nsyst)


struct MetaData
 name::String
 Ld_ePp :: Float64
 Ld_eMp :: Float64
 sqrtS::Float64
end


const ZEUS_MD = MetaData("ZEUS", 141.44, 185.018, 318.0);


"""
    SysPenality(x)

returns gaussian penality for introducing systematics
"""

function SysPenality(x)
    val = (1/sqrt(2*pi))*exp(-(x*x/2.))
    return val
end

"""
    Init_sys()

Reads various systematic errors and feeds them for further use
"""
function Init_sys()
# NP1 429
# NP2 153
# nsyst 8

TM_Elements_ePp  = zeros(NP1,NP2)
TM_Elements_eMp  = zeros(NP1,NP2)
    
 TM_Elements_ePp = get_TM_elements(0,MD_ZEUS);
 TM_Elements_eMp = get_TM_elements(1,MD_ZEUS);


 Tnm_Ee_sys_ePp = zeros(NP1,NP2)
 Tnm_Eehigh_ePp[NP1,NP2]
    
 Tnm_sys_ePp = zeros(NP1,NP2,nsyst)
 Tnm_sys_eMp = zeros(NP1,NP2,nsyst)


    
    TM_Elements_ePp = TM_Elements_ePp / MD_ZEUS.Ld_ePp
    TM_Elements_eMp = TM_Elements_eMp / MD_ZEUS.Ld_eMp
    

 for i in 1:NP1
    for j in 1:NP2

            

     Tnm_Ee_sys_ePp[i,j] =( abs(Tnm_Eehigh_ePp[i,j] - TM_Elements_ePp[i,j])
      +abs(Tnm_Eelow_ePp[i,j] - TM_Elements_ePp[i,j]) )/2.0

     Tnm_sys_ePp[i,j,1] =( abs(Tnm_Eehigh_ePp[i,j] - TM_Elements_ePp[i,j])
      +abs(Tnm_Eelow_ePp[i,j] - TM_Elements_ePp[i,j]) )/2.0

     Tnm_sys_ePp[i,j,2] =( abs(Tnm_Eeconehigh_ePp[i,j] - TM_Elements_ePp[i,j])
      +abs(Tnm_Eeconelow_ePp[i,j] - TM_Elements_ePp[i,j]) )/2.

     Tnm_sys_ePp[i,j,3] =(abs(Tnm_Eereshigh_ePp[i,j] - TM_Elements_ePp[i,j])
      +abs(Tnm_Eereslow_ePp[i,j] - TM_Elements_ePp[i,j]) )/2.


     Tnm_sys_ePp[i,j,4] =( abs(Tnm_Ejhigh_ePp[i,j] - TM_Elements_ePp[i,j])
      +abs(Tnm_Ejlow_ePp[i,j] - TM_Elements_ePp[i,j]))/2.


     Tnm_sys_ePp[i,j,5] =(abs(Tnm_FBcalhigh_ePp[i,j] - TM_Elements_ePp[i,j])
      +abs(Tnm_FBcallow_ePp[i,j] - TM_Elements_ePp[i,j]))/2.

     Tnm_sys_ePp[i,j,6] =(abs(Tnm_FCALxhigh_ePp[i,j] - TM_Elements_ePp[i,j])
      +abs(Tnm_FCALxlow_ePp[i,j] - TM_Elements_ePp[i,j]) )/2.

     Tnm_sys_ePp[i,j,7] =(abs(Tnm_FCALyhigh_ePp[i,j] - TM_Elements_ePp[i,j])
      +abs(Tnm_FCALylow_ePp[i,j] - TM_Elements_ePp[i,j]) )/2.


     Tnm_sys_ePp[i,j,8] =(abs(Tnm_AriMepsUp_ePp[i,j] - TM_Elements_ePp[i,j])
      +abs(Tnm_AriMepsDown_ePp[i,j] - TM_Elements_ePp[i,j]) )/2.




     Tnm_sys_eMp[i,j,1] =(abs(Tnm_Eehigh_eMp[i,j] - TM_Elements_eMp[i,j])
      +abs(Tnm_Eelow_eMp[i,j] - TM_Elements_eMp[i,j]))/2.

     Tnm_sys_eMp[i,j,2] =(abs(Tnm_Eeconehigh_eMp[i,j] - TM_Elements_eMp[i,j])
      +abs(Tnm_Eeconelow_eMp[i,j] - TM_Elements_eMp[i,j]))/2.

     Tnm_sys_eMp[i,j,3] =(abs(Tnm_Eereshigh_eMp[i,j] - TM_Elements_eMp[i,j])
      +abs(Tnm_Eereslow_eMp[i,j] - TM_Elements_eMp[i,j]))/2.


     Tnm_sys_eMp[i,j,4] =(abs(Tnm_Ejhigh_eMp[i,j] - TM_Elements_eMp[i,j])
      +abs(Tnm_Ejlow_eMp[i,j] - TM_Elements_eMp[i,j]) )/2.


     Tnm_sys_eMp[i,j,5] =(abs(Tnm_FBcalhigh_eMp[i,j] - TM_Elements_eMp[i,j])
      +abs(Tnm_FBcallow_eMp[i,j] - TM_Elements_eMp[i,j]) )/2.

     Tnm_sys_eMp[i,j,6] =(abs(Tnm_FCALxhigh_eMp[i,j] - TM_Elements_eMp[i,j])
      +abs(Tnm_FCALxlow_eMp[i,j] - TM_Elements_eMp[i,j]) )/2.

     Tnm_sys_eMp[i,j,7] =(abs(Tnm_FCALyhigh_eMp[i,j] - TM_Elements_eMp[i,j])
      +abs(Tnm_FCALylow_eMp[i,j] - TM_Elements_eMp[i,j]) )/2.


     Tnm_sys_eMp[i,j,8] =(abs(Tnm_AriMepsDown_eMp[i,j] - TM_Elements_eMp[i,j])
      +abs(Tnm_AriMepsUp_eMp[i,j] - TM_Elements_eMp[i,j]) )/2.


    end
 end
    
     println(TM_Elements_ePp[NP1,NP2])
     println(Tnm_Eehigh_ePp[NP1,NP2])
     println(Tnm_Eelow_ePp[NP1,NP2])
     println(Tnm_Ee_sys_ePp[NP1,NP2])
     println(Tnm_sys_ePp[NP1,NP2,1])
    
     Tnm_sys_ePp = Tnm_sys_ePp * MD_ZEUS.Ld_ePp
     Tnm_sys_eMp = Tnm_sys_eMp * MD_ZEUS.Ld_eMp
    
     println(Tnm_Ee_sys_ePp[NP1,NP2])
     println(Tnm_sys_ePp[NP1,NP2,1])
    
end
