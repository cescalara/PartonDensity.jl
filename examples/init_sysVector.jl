#%%

include("./data/ePp_jl/SysTnm_Eehigh_ePp.jl")
include("./data/ePp_jl/SysTnm_Eelow_ePp.jl")
include("./data/ePp_jl/SysTnm_Eeconehigh_ePp.jl")
include("./data/ePp_jl/SysTnm_Eeconelow_ePp.jl")
include("./data/ePp_jl/SysTnm_Eereshigh_ePp.jl")
include("./data/ePp_jl/SysTnm_Eereslow_ePp.jl")
include("./data/ePp_jl/SysTnm_Ejhigh_ePp.jl")
include("./data/ePp_jl/SysTnm_Ejlow_ePp.jl")
include("./data/ePp_jl/SysTnm_FBcalhigh_ePp.jl")
include("./data/ePp_jl/SysTnm_FBcallow_ePp.jl")
include("./data/ePp_jl/SysTnm_FCALxhigh_ePp.jl")
include("./data/ePp_jl/SysTnm_FCALxlow_ePp.jl")
include("./data/ePp_jl/SysTnm_FCALyhigh_ePp.jl")
include("./data/ePp_jl/SysTnm_FCALylow_ePp.jl")
include("./data/ePp_jl/SysTnm_AriMepsUp_ePp.jl")
include("./data/ePp_jl/SysTnm_AriMepsDown_ePp.jl")

include("./data/eMp_jl/SysTnm_Eehigh_eMp.jl")
include("./data/eMp_jl/SysTnm_Eelow_eMp.jl")
include("./data/eMp_jl/SysTnm_Eeconehigh_eMp.jl")
include("./data/eMp_jl/SysTnm_Eeconelow_eMp.jl")
include("./data/eMp_jl/SysTnm_Eereshigh_eMp.jl")
include("./data/eMp_jl/SysTnm_Eereslow_eMp.jl")
include("./data/eMp_jl/SysTnm_Ejhigh_eMp.jl")
include("./data/eMp_jl/SysTnm_Ejlow_eMp.jl")
include("./data/eMp_jl/SysTnm_FBcalhigh_eMp.jl")
include("./data/eMp_jl/SysTnm_FBcallow_eMp.jl")
include("./data/eMp_jl/SysTnm_FCALxhigh_eMp.jl")
include("./data/eMp_jl/SysTnm_FCALxlow_eMp.jl")
include("./data/eMp_jl/SysTnm_FCALyhigh_eMp.jl")
include("./data/eMp_jl/SysTnm_FCALylow_eMp.jl")
include("./data/eMp_jl/SysTnm_AriMepsUp_eMp.jl")
include("./data/eMp_jl/SysTnm_AriMepsDown_eMp.jl")

using PartonDensity   

#include("/home/ritu/learing_julia/data_T.jl")

 
Tnm_Eehigh_ePp = get_TM_elements(0);
Tnm_Eehigh_eMp = get_TM_elements(1);


Tnm_sys_ePp = zeros(429,153,8)
Tnm_sys_eMp = zeros(429,153,8)
Tnm_Ee_sys_ePp = zeros(429,153)
Tnm_Eehigh_ePp[426,153]

for i in 1:429
    for j in 1:153


     Tnm_Ee_sys_ePp[i,j] =abs(Tnm_Eehigh_ePp[i,j] - TM_Elements_ePp[i,j])/2.
      +abs(Tnm_Eelow_ePp[i,j] - TM_Elements_ePp[i,j])/2.

     Tnm_sys_ePp[i,j,1] =abs(Tnm_Eehigh_ePp[i,j] - TM_Elements_ePp[i,j])/2.
      +abs(Tnm_Eelow_ePp[i,j] - TM_Elements_ePp[i,j])/2.

     Tnm_sys_ePp[i,j,2] =abs(Tnm_Eeconehigh_ePp[i,j] - TM_Elements_ePp[i,j])/2.
      +abs(Tnm_Eeconelow_ePp[i,j] - TM_Elements_ePp[i,j])/2.

     Tnm_sys_ePp[i,j,3] =abs(Tnm_Eereshigh_ePp[i,j] - TM_Elements_ePp[i,j])/2.
      +abs(Tnm_Eereslow_ePp[i,j] - TM_Elements_ePp[i,j])/2.


     Tnm_sys_ePp[i,j,4] =abs(Tnm_Ejhigh_ePp[i,j] - TM_Elements_ePp[i,j])/2.
      +abs(Tnm_Ejlow_ePp[i,j] - TM_Elements_ePp[i,j])/2.


     Tnm_sys_ePp[i,j,5] =abs(Tnm_FBcalhigh_ePp[i,j] - TM_Elements_ePp[i,j])/2.
      +abs(Tnm_FBcallow_ePp[i,j] - TM_Elements_ePp[i,j])/2.

     Tnm_sys_ePp[i,j,6] =abs(Tnm_FCALxhigh_ePp[i,j] - TM_Elements_ePp[i,j])/2.
      +abs(Tnm_FCALxlow_ePp[i,j] - TM_Elements_ePp[i,j])/2.

     Tnm_sys_ePp[i,j,7] =abs(Tnm_FCALyhigh_ePp[i,j] - TM_Elements_ePp[i,j])/2.
      +abs(Tnm_FCALylow_ePp[i,j] - TM_Elements_ePp[i,j])/2.


     Tnm_sys_ePp[i,j,8] =abs(Tnm_AriMepsUp_ePp[i,j] - TM_Elements_ePp[i,j])/2.
      +abs(Tnm_AriMepsDown_ePp[i,j] - TM_Elements_ePp[i,j])/2.




     Tnm_sys_eMp[i,j,1] =abs(Tnm_Eehigh_eMp[i,j] - TM_Elements_eMp[i,j])/2.
      +abs(Tnm_Eelow_eMp[i,j] - TM_Elements_eMp[i,j])/2.

     Tnm_sys_eMp[i,j,2] =abs(Tnm_Eeconehigh_eMp[i,j] - TM_Elements_eMp[i,j])/2.
      +abs(Tnm_Eeconelow_eMp[i,j] - TM_Elements_eMp[i,j])/2.

     Tnm_sys_eMp[i,j,3] =abs(Tnm_Eereshigh_eMp[i,j] - TM_Elements_eMp[i,j])/2.
      +abs(Tnm_Eereslow_eMp[i,j] - TM_Elements_eMp[i,j])/2.


     Tnm_sys_eMp[i,j,4] =abs(Tnm_Ejhigh_eMp[i,j] - TM_Elements_eMp[i,j])/2.
      +abs(Tnm_Ejlow_eMp[i,j] - TM_Elements_eMp[i,j])/2.


     Tnm_sys_eMp[i,j,5] =abs(Tnm_FBcalhigh_eMp[i,j] - TM_Elements_eMp[i,j])/2.
      +abs(Tnm_FBcallow_eMp[i,j] - TM_Elements_eMp[i,j])/2.

     Tnm_sys_eMp[i,j,6] =abs(Tnm_FCALxhigh_eMp[i,j] - TM_Elements_eMp[i,j])/2.
      +abs(Tnm_FCALxlow_eMp[i,j] - TM_Elements_eMp[i,j])/2.

     Tnm_sys_eMp[i,j,7] =abs(Tnm_FCALyhigh_eMp[i,j] - TM_Elements_eMp[i,j])/2.
      +abs(Tnm_FCALylow_eMp[i,j] - TM_Elements_eMp[i,j])/2.


     Tnm_sys_eMp[i,j,8] =abs(Tnm_AriMepsDown_eMp[i,j] - TM_Elements_eMp[i,j])/2.
      +abs(Tnm_AriMepsUp_eMp[i,j] - TM_Elements_eMp[i,j])/2.


    end
end

# check reading the sys variation at one point

Tnm_Ee_sys_ePp[429,153]
Tnm_sys_ePp[429,153,1]

