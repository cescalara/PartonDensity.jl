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


include("../data/zeus_transfer_matrix.jl")

export get_pred_N
export get_data_events


function get_data_events(eMPp)
    if eMPp == 0
        return Data_Events_ePp
    else
        return Data_Events_eMp
    end
end


function get_K_elements(eMPp)
    if eMPp == 0
        return Kii_ePp
    else
        return Kii_eMp
    end
end


function get_L_data(eMPp)
    if eMPp == 0
        return Ld_ePp
    else
        return Ld_eMp
    end
end

function get_TM_elements(eMPp)
    lumi_data = get_L_data(eMPp)
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

function get_pred_N(Integ_xsec, eMPp)

    TM = get_TM_elements(eMPp);
    K = get_K_elements(eMPp)

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
