# # ZEUS interface
#
# Here is a short demonstration of how to access the ZEUS transfer matrix and bin interface.

using PartonDensity, CSV, DelimitedFiles

"""
    get_bin_info(n, quiet)

Get the bin egdes of the ZEUS 
detector space for a given 
bin number, `n`.
"""
function get_bin_info(n::Integer; quiet::Bool = false)

    if n < 1 || n > 153
        @error "Bin number n should be [1, 153]"
    end
    if !quiet
        @info "ZEUS detector bin" n m_BinQ2low[n] m_BinQ2high[n] m_Binxlow[n] m_Binxhigh[n]
    end
    return ([m_BinQ2low[n], m_BinQ2high[n]], [m_Binxlow[n], m_Binxhigh[n]])
end

function get_pred_N(Integ_xsec, eMPp, md::MetaData)
    if eMPp == 0
    TM = m_TM_elements_ePp
    K = m_K_elements_ePp
    else
    TM = m_TM_elements_eMp
    K = m_K_elements_eMp)    
    end
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


# ## Transfer matrix

eMPp = 1 # e+/e- switch 0/1

# Read in an example integrated cross section
numbers_from_file = readdlm("data/EXAMPLE_1/HERAPDF20_NNLO_EIG_ePp.txt") 

# List of integrated cross section values in 429 bins 
integ_xsec = numbers_from_file[:,3] 

# Corresponding list of expected event numbers
prediction = get_pred_N(integ_xsec, eMPp, MD_ZEUS_I1787035); 

integ_xsec[153]

prediction[151]

# ## Bins in detector space
#
# The transfer matrix projects into a list of 153 bins with irregular x and Q2 spacing.
# The bin edges are provided in `data/ZEUS_I1787035/ZEUS_I1787035.jl` but can also be quickly and
# easily accessed with the helper function `get_bin_info` as shown below.

get_bin_info(10)
