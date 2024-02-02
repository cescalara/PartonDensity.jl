# # ZEUS interface
#
# Here is a short demonstration of how to access the ZEUS transfer matrix and bin interface.

using PartonDensity, CSV, DelimitedFiles

zeus_include_path = string(chop(pathof(PartonDensity), tail=20), "data/ZEUS_I1787035/ZEUS_I1787035.jl")
sample_data_include_path = string(chop(pathof(PartonDensity), tail=20), "data/EXAMPLE_1/HERAPDF20_NNLO_EIG_ePp.txt")
MD_ZEUS_I1787035=include(zeus_include_path)

# Read in an example integrated cross section
numbers_from_file = readdlm(sample_data_include_path)

# List of integrated cross section values in 429 bins 
integ_xsec_ePp = numbers_from_file[:, 3]
integ_xsec_eMp = numbers_from_file[:, 3]

# Corresponding list of expected event numbers
prediction_ePp, prediction_eMp = MD_ZEUS_I1787035.f_cross_section_to_counts(MD_ZEUS_I1787035.Ld_ePp, MD_ZEUS_I1787035.Ld_eMp, integ_xsec_ePp, integ_xsec_eMp)

integ_xsec_ePp[153]

prediction_ePp[151]

