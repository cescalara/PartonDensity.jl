#%%
module MyModule

"""
// ****************************************
// R. Aggarwal
// Tranfer Matrix Analysis for the high-x data
// should be used as N = T K S
// K - radiative corrections
// S - Born level integrated sigma .. to be read externally
//       (Units pb-1)
// T- Transfer Matrix (has data luminosity taken into account)
//******************************************* 
"""


include("data_T.jl")

export Get_Pred_N
export Get_data_events



function Get_data_events(eMPp)
    if eMPp == 0
        return Data_Events_ePp
    else
        return Data_Events_eMp
    end
end


function Get_K_Elements(eMPp)
    if eMPp == 0
        return Kii_ePp
    else
        return Kii_eMp
    end
end


function Get_L_data(eMPp)
    if eMPp == 0
        return Ld_ePp
    else
        return Ld_eMp
    end
end

function Get_TM_Elements(eMPp)
    lumi_data = Get_L_data(eMPp)
    TM_elements = TM_Elements_eMp
    if eMPp == 0
        TM_elements = TM_Elements_ePp
    end
    for i in 1:429
        for j in 1:153
         TM_elements[i,j] =TM_elements[i,j] * lumi_data
        end
    end
    return TM_elements
end

function Get_Pred_N(Integ_xsec, eMPp)

    TM = Get_TM_Elements(eMPp);
    K=Get_K_Elements(eMPp)

    xsec_pred = zeros(153)
    for j in 1:153
        temp=0.
        for i in 1:429

        temp = temp +(TM[i,j]./K[i]) * Integ_xsec[i];
        end
        xsec_pred[j] = temp
    end

    return xsec_pred;

end

end #end of module

#%%