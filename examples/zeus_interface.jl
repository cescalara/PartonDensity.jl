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


# Read in an example integrated cross section
numbers_from_file = readdlm("data/EXAMPLE_1/HERAPDF20_NNLO_EIG_ePp.txt") 

# List of integrated cross section values in 429 bins 
integ_xsec_ePp = numbers_from_file[:,3] 
integ_xsec_eMp = numbers_from_file[:,3] 

# Corresponding list of expected event numbers
prediction_ePp, prediction_eMp = f_cross_section_to_counts(integ_xsec_ePp,integ_xsec_eMp)

integ_xsec_ePp[153]

prediction_ePp[151]

# ## Bins in detector space
#
# The transfer matrix projects into a list of 153 bins with irregular x and Q2 spacing.
# The bin edges are provided in `data/ZEUS_I1787035/ZEUS_I1787035.jl` but can also be quickly and
# easily accessed with the helper function `get_bin_info` as shown below.

get_bin_info(10)
