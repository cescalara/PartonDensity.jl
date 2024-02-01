# Introduction

A Bayesian approach to parton density extraction. 

Hadrons, such as protons and neutrons, are made up of quarks held together by the strong force. At high energy scales, the valence quarks that define these hadrons exist in a sea of virtual quarks and gluons. The parton distribution functions (PDFs) describe this structure and are of fundamental importance to our understanding of quantum chromodynamics (QCD), as well as its application to LHC physics and the development of cosmic ray air showers in the Earth's atmosphere. PDFs can be extracted from accelerator measurements in which hadrons are probed through collisions with electrons. A limitation of existing approaches to analysing this data is the reliance on the chi-square statistic and the coupled assumption of Normal-distributed observations. We are working on a new statistical method for PDF extraction, which overcomes this limitation by forward modelling the problem from an input PDF to the expected number of events in a detector. This model will then be fit using Markov Chain Monte Carlo to enable inference of the PDF parameters. Our project builds on the QCDNUM software for fast QCD evolution and the Bayesian Analysis Toolkit for inference. We initially focus on the "high-x" regime and data from the [ZEUS experiment](https://particle-physics.desy.de/research/previous_desy_experiments/zeus/index_eng.html), where the chi-square method cannot be used due to low event numbers.

This package uses [QCDNUM.jl](https://github.com/cescalara/QCDNUM.jl) for fast PDF evolution and cross-section calculation and [BAT.jl](https://github.com/bat/BAT.jl) for Bayesian inference.

PartonDensity.jl has been used in the following publications:
* [arXiv:2209.06571](https://arxiv.org/abs/2209.06571)
* [arXiv:2401.17729](https://arxiv.org/abs/2401.17729)

## Installation

To install PartonDensity.jl, start Julia and run

```julia
julia> using Pkg
julia> Pkg.add(url="https://github.com/cescalara/PartonDensity.jl.git")
```

## Usage

Check out the examples listed in these docs! The scripts can be found in the [examples](https://github.com/cescalara/PartonDensity.jl/tree/main/examples) directory of the GitHub repository. 

To run these docs, you can follow the steps used in the GitHub workflow. In particular, we need to remember to use the latest versions of `BAT.jl` from GitHub. 

```julia
using Pkg
Pkg.add(url="https://github.com/bat/BAT.jl.git")
```

If you want to convert the example scripts to notebooks, use `Literate.jl`.

## Development

Below are the installation instruction for those who wish to contribute to the code.

- Clone the github repository, e.g. via the command line:
```
git clone  https://github.com/cescalara/PartonDensity.jl.git
```

- Enter the directory and start Julia interpreter
```
cd PartonDensity.jl
julia
```

-  Open the Julia package management environment pressing `]`.

```
julia> ]
```

 - Execute 
```
pkg> generate PartonDensity
...... 
pkg>  . dev
```
 - Exit the package manager using backspace or pressing `Ctrl+C`

