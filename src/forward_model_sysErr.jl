using HDF5

export forward_model_sysErr, forward_model_init_sysErr, SysPenality
export Tnm_sys_ePp, Tnm_sys_eMp

include("../examples/data/ePp_jl/SysTnm_Eehigh_ePp.jl")
include("../examples/data/ePp_jl/SysTnm_Eelow_ePp.jl")
include("../examples/data/ePp_jl/SysTnm_Eeconehigh_ePp.jl")
include("../examples/data/ePp_jl/SysTnm_Eeconelow_ePp.jl")
include("../examples/data/ePp_jl/SysTnm_Eereshigh_ePp.jl")
include("../examples/data/ePp_jl/SysTnm_Eereslow_ePp.jl")
include("../examples/data/ePp_jl/SysTnm_Ejhigh_ePp.jl")
include("../examples/data/ePp_jl/SysTnm_Ejlow_ePp.jl")
include("../examples/data/ePp_jl/SysTnm_FBcalhigh_ePp.jl")
include("../examples/data/ePp_jl/SysTnm_FBcallow_ePp.jl")
include("../examples/data/ePp_jl/SysTnm_FCALxhigh_ePp.jl")
include("../examples/data/ePp_jl/SysTnm_FCALxlow_ePp.jl")
include("../examples/data/ePp_jl/SysTnm_FCALyhigh_ePp.jl")
include("../examples/data/ePp_jl/SysTnm_FCALylow_ePp.jl")
include("../examples/data/ePp_jl/SysTnm_AriMepsUp_ePp.jl")
include("../examples/data/ePp_jl/SysTnm_AriMepsDown_ePp.jl")

include("../examples/data/eMp_jl/SysTnm_Eehigh_eMp.jl")
include("../examples/data/eMp_jl/SysTnm_Eelow_eMp.jl")
include("../examples/data/eMp_jl/SysTnm_Eeconehigh_eMp.jl")
include("../examples/data/eMp_jl/SysTnm_Eeconelow_eMp.jl")
include("../examples/data/eMp_jl/SysTnm_Eereshigh_eMp.jl")
include("../examples/data/eMp_jl/SysTnm_Eereslow_eMp.jl")
include("../examples/data/eMp_jl/SysTnm_Ejhigh_eMp.jl")
include("../examples/data/eMp_jl/SysTnm_Ejlow_eMp.jl")
include("../examples/data/eMp_jl/SysTnm_FBcalhigh_eMp.jl")
include("../examples/data/eMp_jl/SysTnm_FBcallow_eMp.jl")
include("../examples/data/eMp_jl/SysTnm_FCALxhigh_eMp.jl")
include("../examples/data/eMp_jl/SysTnm_FCALxlow_eMp.jl")
include("../examples/data/eMp_jl/SysTnm_FCALyhigh_eMp.jl")
include("../examples/data/eMp_jl/SysTnm_FCALylow_eMp.jl")
include("../examples/data/eMp_jl/SysTnm_AriMepsUp_eMp.jl")
include("../examples/data/eMp_jl/SysTnm_AriMepsDown_eMp.jl")
include("../examples/data/zeus_bin_edges.jl")


UInt64::NP1 = size(BinQ2low)
UInt64::NP2 = size(Binxlow)
UInt64::nsyst = 8


Tnm_sys_ePp = zeros(NP1,NP1,nsyst)
Tnm_sys_eMp = zeros(NP1,NP2,nsyst)


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
    
 TM_Elements_ePp = get_TM_elements(0);
 TM_Elements_eMp = get_TM_elements(1);


 Tnm_Ee_sys_ePp = zeros(NP1,NP2)
 Tnm_Eehigh_ePp[NP1,NP2]
    
 Tnm_sys_ePp = zeros(NP1,NP2,nsyst)
 Tnm_sys_eMp = zeros(NP1,NP2,nsyst)


    
    TM_Elements_ePp = TM_Elements_ePp / get_L_data(0)
    TM_Elements_eMp = TM_Elements_eMp / get_L_data(1)
    

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
    
     Tnm_sys_ePp = Tnm_sys_ePp * get_L_data(0)
     Tnm_sys_eMp = Tnm_sys_eMp * get_L_data(1)
    
     println(Tnm_Ee_sys_ePp[NP1,NP2])
     println(Tnm_sys_ePp[NP1,NP2,1])
    
end

"""
    forward_model_init(qcdnum_grid, qcdnum_params)

Initialise forward model with systematic Error.
Initialises QCDNUM and builds weight tables to 
save time in subsequent iterations. Must be called prior to first instance 
of `forward_model()`.
"""
function forward_model_init_sysErr(qcdnum_grid::QCDNUMGrid, qcdnum_params::QCDNUMParameters,
                            splint_params::SPLINTParameters)

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
    # TODO: Use saved weights file once QCDNUM fixed
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
    Init_sys()
   
end

