# # ZEUS interface
#
# Here is a short demonstration of how to access the ZEUS transfer matrix and bin interface.

using PartonDensity, CSV, DelimitedFiles

# ## Transfer matrix

eMPp = 1 # e+/e- switch 0/1

# Read in an example integrated cross section
numbers_from_file = readdlm("../data/EXAMPLE_1/HERAPDF20_NNLO_EIG_ePp.txt") 

# List of integrated cross section values in 429 bins 
integ_xsec = numbers_from_file[:,3] 

# Corresponding list of expected event numbers
prediction = get_pred_N(integ_xsec, eMPp, MD_ZEUS); 

integ_xsec[153]

prediction[151]

# ## Bins in detector space
#
# The transfer matrix projects into a list of 153 bins with irregular x and Q2 spacing.
# The bin edges are provided in `data/ZEUS_I1787035/ZEUS_I1787035.jl` but can also be quickly and
# easily accessed with the helper function `get_bin_info` as shown below.

get_bin_info(10)
