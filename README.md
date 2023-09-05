# PartonDensity

[![](https://img.shields.io/badge/docs-dev-blue.svg)](https://francescacapel.com/PartonDensity.jl/dev/) 
![CI](https://github.com/cescalara/PartonDensity.jl/actions/workflows/Testsv2.yml/badge.svg)
[![codecov](https://codecov.io/gh/cescalara/PartonDensity.jl/branch/main/graph/badge.svg?token=Q5PLHWXV3I)](https://codecov.io/gh/cescalara/PartonDensity.jl)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

A Bayesian approach to parton density extraction.

For more information please see the [documentation](https://francescacapel.com/PartonDensity.jl/) and `examples/` scripts.

## Quick start

To install PartonDensity.jl and its dependencies, start Julia and run

```julia
julia> using Pkg
julia> Pkg.add(url="https://github.com/cescalara/PartonDensity.jl.git")
```
To install a specific version (git tag) of PartonDensity.jl, e.g. "v0.1.0", run
```
julia> Pkg.add(url="https://github.com/cescalara/PartonDensity.jl.git", rev="v0.1.0")
```

To run the examples scripts, you can follow the steps used in the GitHub workflow. In particular, we need to remember to use the latest versions of ``BAT.jl`` from GitHub. 

```julia
using Pkg
Pkg.add(url="https://github.com/bat/BAT.jl.git")
```
