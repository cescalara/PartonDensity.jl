---
jupyter:
  jupytext:
    text_representation:
      extension: .md
      format_name: markdown
      format_version: '1.3'
      jupytext_version: 1.11.0
  kernelspec:
    display_name: Julia 1.6.1
    language: julia
    name: julia-1.6
---

```julia
using PartonDensity, CSV, DelimitedFiles
```

```julia
eMPp = 1 # e+/e- switch 0/1
numbers_from_file = readdlm("data/HERAPDF20_NNLO_EIG_ePp.txt") 

# List of integrated cross section values in 429 bins 
integ_xsec = numbers_from_file[:,3] 

# Corresponding list of expected event numbers
prediction = get_pred_N(integ_xsec, eMPp); 
```

```julia
integ_xsec[153]
```

```julia
prediction[151]
```

```julia

```
