# # Transfer matrix
#
# Here is a short demonstration of how to access the ZEUS transfer matrix interface.

using PartonDensity, CSV, DelimitedFiles

eMPp = 1 # e+/e- switch 0/1

# Read in an example integrated cross section
numbers_from_file = readdlm("data/HERAPDF20_NNLO_EIG_ePp.txt") 

# List of integrated cross section values in 429 bins 
integ_xsec = numbers_from_file[:,3] 

# Corresponding list of expected event numbers
prediction = get_pred_N(integ_xsec, eMPp); 

integ_xsec[153]

prediction[151]
