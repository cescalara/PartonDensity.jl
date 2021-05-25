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
eMPp = 1
numbers = readdlm("data/HERAPDF20_NNLO_EIG_ePp.txt")
Intig = numbers[:,3]
prediction = Get_Pred_N(Intig,eMPp);
```

```julia
Intig[153]
```

```julia
prediction[151]
```

```julia

```
