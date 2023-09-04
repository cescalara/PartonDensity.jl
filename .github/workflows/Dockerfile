FROM fedora:38
ENV LC_ALL=C
ENV JULIA_PKGDIR=/opt/julia
ENV JULIA_DEPOT_PATH=/opt/myjulia
RUN  yum -y  install  dnf-plugins-core \
                      bc make cmake binutils git wget diffutils file sed gawk grep which autoconf automake libtool \
                      gcc-gfortran gcc-c++ clang bzip2   openssl-devel openssl \
                      julia python  && yum -y clean all && \
julia  -e 'import Pkg;                                 \
Pkg.add("SpecialFunctions");                           \
Pkg.add(url="https://github.com/bat/BAT.jl.git");      \
Pkg.add("Colors");                                     \
Pkg.add("Plots");                                      \
Pkg.add("HDF5"); \
Pkg.add("PackageCompiler"); \
Pkg.add("ArgParse"); \
Pkg.add("LaTeXStrings"); \
Pkg.add("DensityInterface"); \
Pkg.add("QCDNUM"); \
Pkg.add("Plots"); \
Pkg.add("Random"); \
Pkg.add("Distributions"); \
Pkg.add("ValueShapes"); \
Pkg.add("ParallelProcessingTools"); \
Pkg.add("StatsBase"); \
Pkg.add("LinearAlgebra"); \
Pkg.add("SpecialFunctions"); \
Pkg.add("Printf"); \
Pkg.add("DelimitedFiles"); \
Pkg.add("LaTeXStrings"); \
Pkg.add("HypothesisTests"); \
Pkg.add("Statistics"); \
Pkg.add("Measures"); \
Pkg.add("ArgParse");'
