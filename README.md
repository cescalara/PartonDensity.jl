# PartonDensity

[![](https://img.shields.io/badge/docs-dev-blue.svg)](https://francescacapel.com/PartonDensity.jl/dev/) 
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

A Bayesian approach to parton density extraction.

For more information please see the [documentation](https://francescacapel.com/PartonDensity.jl/) and `examples/` scripts.

## Quick start

To install PartonDensity.jl and its dependencies, start Julia and run

```julia
julia> using Pkg
julia> Pkg.add(url="https://github.com/cescalara/PartonDensity.jl.git")
```

To run the examples scripts, you can follow the steps used in the GitHub workflow. In particular, we need to remember to use the latest versions of ``QCDNUM.jl`` and ``BAT.jl`` from GitHub. 

```julia
using Pkg
Pkg.add(url="https://github.com/bat/BAT.jl.git")
Pkg.add(url="https://github.com/cescalara/QCDNUM.jl.git")
```
