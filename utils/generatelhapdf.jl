#!/usr/bin/julia
using BAT, DensityInterface
using PartonDensity
using QCDNUM
using Plots, Colors , Random, Distributions, ValueShapes, ParallelProcessingTools
using StatsBase, LinearAlgebra
using Printf
using ArgParse
pdf_params = DirichletPDFParams(K_u=3.1, K_d=3.7, λ_g1=0.7, λ_g2=-0.5, K_g=5.0,λ_q=-0.5, K_q=6.0,θ=[0.22, 0.10, 0.24, 0.24, 0.10,0.05, 0.01, 0.005, 0.0005])

function func(ipdf, x)::Float64
#pdf_params = DirichletPDFParams(K_u=3.7, K_d=3.7, λ_g1=0.5, λ_g2=-0.5, K_g=5.0,λ_q=-0.5, K_q=6.0,θ=[0.22, 0.10, 0.24, 0.24, 0.10,0.05, 0.01, 0.005, 0.0005])
r::Float64 = get_input_pdf_func(pdf_params)(ipdf, x)
return r
end

function parse_commandline()
    s = ArgParseSettings()

    @add_arg_table s begin
        "--seed", "-s"
            help = "Seed"
            arg_type = Int
            default = 42
        "--parametrisation", "-p"
            help = "Parametrisation -- Dirichlet or Valence"
            arg_type = String
            default = "Dirichlet"
        "--fitresults", "-f"
            help = "Input fitresults -- file in the pseudodata directory w/o the extension"
            arg_type = String
            default = "fit-Dirichlet-0-42-data"    
            
            
    end

    return parse_args(s)
end


function main()
parsed_args = parse_commandline()
println("Parsed args:")
for (arg,val) in parsed_args
  println("  $arg  =>  $val")
end
seed=parsed_args["seed"]
println(seed)
seedtxt=string(seed)
rng = MersenneTwister(seed)

qcdnum_grid = QCDNUM.GridParams(x_min=[1.0e-3, 1.0e-1, 5.0e-1], x_weights=[1, 2, 2], nx=100, qq_bounds=[1.0e2, 3.0e4], qq_weights=[1.0, 1.0], nq=50, spline_interp=3)
qcdnum_params = QCDNUM.EvolutionParams(order=2, α_S=0.118, q0=100.0, grid_params=qcdnum_grid,n_fixed_flav=5, iqc=1, iqb=1, iqt=1, weight_type=1);
    
splint_params = QCDNUM.SPLINTParams();
quark_coeffs = QuarkCoefficients();
forward_model_init(qcdnum_params, splint_params)

func_c = @cfunction(func, Float64, (Ref{Int32}, Ref{Float64}))


xmin = Float64.([1.0e-4])
iwt = Int32.([1])
ng = 1
nxin = 100
iosp = 3
nx = 10

qq = Float64.([2e0, 1e4])
wt = Float64.([1e0, 1e0])
nq = 1
nqin = 60
ngq = 2
itype = 1

as0 = 0.364
r20 = 2.0

q2c = 3.0
q2b = 25.0
q0 = 2.0
#q0 = 100.0
iqt = 999

def = Float64.([0., 0., 0., 0., 0.,-1., 0., 1., 0., 0., 0., 0., 0.,     
                0., 0., 0., 0.,-1., 0., 0., 0., 1., 0., 0., 0., 0.,      
                0., 0., 0.,-1., 0., 0., 0., 0., 0., 1., 0., 0., 0.,      
                0., 0., 0., 0., 0., 1., 0., 0., 0., 0., 0., 0., 0.,      
                0., 0., 0., 0., 1., 0., 0., 0., 0., 0., 0., 0., 0.,      
                0., 0., 0., 1., 0., 0., 0., 0., 0., 0., 0., 0., 0.,      
                0., 0.,-1., 0., 0., 0., 0., 0., 0., 0., 1., 0., 0.,      
                0., 0., 1., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,      
                0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,      
                0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,      
                0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,      
                0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.]);

nfin = 0
x = 1.0e-3
q = 1.0e3
pdf = Array{Float64}(undef, 13)

g = QCDNUM.load_params(string("fitresults/", parsed_args["fitresults"], ".h5"))



QCDNUM.qcinit(-6, " ")
nx = QCDNUM.gxmake(xmin, iwt, ng, nxin, iosp)
nq = QCDNUM.gqmake(qq, wt, ngq, nqin)
nw = QCDNUM.fillwt(itype)
QCDNUM.setord(3)
QCDNUM.setalf(as0, r20)
iqc = QCDNUM.iqfrmq(q2c)
iqb = QCDNUM.iqfrmq(q2b)
iqt = QCDNUM.iqfrmq(1e11)
QCDNUM.setcbt(0, iqc, iqb, 0)
iq0 = QCDNUM.iqfrmq(q0)
#println("OK1")
QCDNUM.evolfg(itype, func_c, def, iq0)

allx = Float64.([
1.113780E-03,1.360370E-03,1.661560E-03,2.029430E-03,2.478750E-03,3.027550E-03,
3.697860E-03,4.516580E-03,5.516560E-03,6.737950E-03,8.229750E-03,1.005180E-02,
1.116455E-02,1.227730E-02,1.363645E-02,1.499560E-02,1.665560E-02,1.831560E-02,
2.034320E-02,2.237080E-02,2.484725E-02,2.732370E-02,3.034850E-02,3.337330E-02,
3.706775E-02,4.076220E-02,4.527465E-02,4.978710E-02,5.529860E-02,6.081010E-02,
6.754185E-02,7.427360E-02,8.249580E-02,9.071800E-02,1.007605E-01,1.108030E-01,
1.230690E-01,1.353350E-01,1.400000E-01,1.446650E-01,1.496515E-01,1.546380E-01,
1.599685E-01,1.652990E-01,1.709965E-01,1.766940E-01,1.827850E-01,1.888760E-01,
1.953865E-01,2.018970E-01,2.088560E-01,2.158150E-01,2.232540E-01,2.306930E-01,
2.386450E-01,2.465970E-01,2.550970E-01,2.635970E-01,2.726830E-01,2.817690E-01,
2.914815E-01,3.011940E-01,3.115760E-01,3.219580E-01,3.330560E-01,3.441540E-01,
3.560165E-01,3.678790E-01,3.805600E-01,3.932410E-01,4.067955E-01,4.203500E-01,
4.348395E-01,4.493290E-01,4.648170E-01,4.803050E-01,4.968610E-01,5.134170E-01,
5.311145E-01,5.488120E-01,5.549780E-01,5.611440E-01,5.674485E-01,5.737530E-01,
5.801995E-01,5.866460E-01,5.932375E-01,5.998290E-01,6.065680E-01,6.133070E-01,
6.201980E-01,6.270890E-01,6.341345E-01,6.411800E-01,6.483840E-01,6.555880E-01,
6.629540E-01,6.703200E-01,6.778515E-01,6.853830E-01,6.930835E-01,7.007840E-01,
7.086575E-01,7.165310E-01,7.245815E-01,7.326320E-01,7.408635E-01,7.490950E-01,
7.575115E-01,7.659280E-01,7.745335E-01,7.831390E-01,7.919380E-01,8.007370E-01,
8.097340E-01,8.187310E-01,8.217745E-01,8.248180E-01,8.278840E-01,8.309500E-01,
8.340390E-01,8.371280E-01,8.402400E-01,8.433520E-01,8.464875E-01,8.496230E-01,
8.527815E-01,8.559400E-01,8.591215E-01,8.623030E-01,8.655090E-01,8.687150E-01,
8.719440E-01,8.751730E-01,8.784265E-01,8.816800E-01,8.849575E-01,8.882350E-01,
8.915370E-01,8.948390E-01,8.981655E-01,9.014920E-01,9.048435E-01,9.081950E-01,
9.115710E-01,9.149470E-01,9.183485E-01,9.217500E-01,9.251765E-01,9.286030E-01,
9.320550E-01,9.355070E-01,9.389845E-01,9.424620E-01,9.459660E-01,9.494700E-01,
9.529995E-01,9.565290E-01,9.600845E-01,9.636400E-01,9.672225E-01,9.708050E-01,
9.744140E-01,9.780230E-01,9.816585E-01,9.852940E-01,9.889570E-01,9.926200E-01,
9.963100E-01,1.000000E+00])
#allq = Float64.([100.,200,300,400,500,1000,5000,10000])
allq = Float64.([
1.096570E+01,1.388560E+01,1.779290E+01,2.308550E+01,3.034710E+01,4.044480E+01,
5.468640E+01,7.507240E+01,1.047120E+02,1.485170E+02,2.143800E+02,3.152120E+02,
4.725370E+02,7.229460E+02,1.129950E+03,1.806160E+03,2.955930E+03])
allalpha = Float64.([
1.096570E+01,1.388560E+01,1.779290E+01,2.308550E+01,3.034710E+01,4.044480E+01,
5.468640E+01,7.507240E+01,1.047120E+02,1.485170E+02,2.143800E+02,3.152120E+02,
4.725370E+02,7.229460E+02,1.129950E+03,1.806160E+03,2.955930E+03])
for (i, qq) in enumerate(allq)
  allalpha[i] =QCDNUM.asfunc(qq*qq)[1]
end
    Ns = 50

open("CABCHSV2023nnlo/CABCHSV2023nnlo.info", "w") do f
  write(f,"
SetDesc: \"CABCHSV PDF fixed 5-flavour fits\"
SetIndex: 83300
Authors: Francesca Capel, Ritu Aggarwal, Michiel Botje, Allen Caldwell, Richard Hildebrandt, Oliver Schulz, Andrii Verbytskyi
Reference: arXiv:2309.xxxx
Format: lhagrid1
DataVersion: 1
NumMembers: ")
write(f,string(Ns+1))
write(f,"\n")
write(f,"Particle: 2212
Flavors: [-5,-4, -3, -2, -1, 1, 2, 3, 4, 5, 21]
OrderQCD: 2
ForcePositive: 1
FlavorScheme: fixed
NumFlavors: 5
XMin: 0.001
XMax: 1
QMin: 100
QMax: 30000
MZ: 91.1876
MUp: 0
MDown: 0
MStrange: 0
MCharm: 1.3
MBottom: 4.75
MTop: 172
AlphaS_MZ: 0.118
AlphaS_OrderQCD: 2
AlphaS_Type: ipol
AlphaS_Qs: ")
write(f,string(allq))
write(f,"\n")
write(f,"AlphaS_Vals: ")
write(f,string(allalpha) )
write(f,"\n")
write(f,"AlphaS_Lambda4: 0.326
AlphaS_Lambda5: 0.226
\n")

end


    samples_data = bat_read(string("fitresults/", parsed_args["fitresults"], ".h5")).result;
    sub_samples_all = BAT.bat_sample(samples_data, BAT.OrderedResampling(nsamples=100)).result;
    samples_mode2 = mode(sub_samples_all);
    sub_samples = BAT.bat_sample(samples_data, BAT.OrderedResampling(nsamples=Ns)).result;
    pdf_params_s = DirichletPDFParams(K_u=samples_mode2.K_u, K_d=samples_mode2.K_d, K_q=samples_mode2.K_q,
                                      λ_g1=samples_mode2.λ_g1, 
                                      λ_g2=samples_mode2.λ_g2,
                                      K_g=samples_mode2.K_g, 
                                      λ_q=samples_mode2.λ_q, 
                                      θ=Vector(samples_mode2.θ))
    allparams = [ pdf_params_s ]
    NN=0
    for s in eachindex(sub_samples)
      push!(allparams,DirichletPDFParams(K_u=sub_samples.v.K_u[s], K_d=sub_samples.v.K_d[s], K_q=sub_samples.v.K_q[s],
                                      λ_g1=sub_samples.v.λ_g1[s], 
                                      λ_g2=sub_samples.v.λ_g2[s],
                                      K_g=sub_samples.v.K_g[s], 
                                      λ_q=sub_samples.v.λ_q[s], 
                                      θ=Vector(sub_samples.v.θ[s]))
      )
    end
    for s in  allparams
      global pdf_params=s
      println(pdf_params)
      QCDNUM.evolfg(itype, func_c, def, iq0)
      open(string("CABCHSV2023nnlo/CABCHSV2023nnlo_", lpad(string(NN),4,"0"),".dat"), "w") do f
      if (NN==0)
        write(f,"PdfType: central\nFormat: lhagrid1\n---\n")
      else
        write(f,"PdfType: error\nFormat: lhagrid1\n---\n")
      end
      NN=NN+1
      write(f,["$v " for v in allx]...) 
      write(f,"\n")
      write(f,["$v " for v in allq]...) 
      write(f,"\n")
      write(f,"-5 -4 -3 -2 -1 1 2 3 4 5 21")
      write(f,"\n")
      for xx in allx
        for qq in allq
          #@printf("%f\n",qq)
          pdf = QCDNUM.allfxq(itype, xx, qq, 0, 1)
          @printf(f, "%8.8f %8.8f %8.8f %8.8f %8.8f %8.8f %8.8f %8.8f %8.8f %8.8f %8.8f\n",pdf[1],pdf[2],pdf[3],pdf[4],pdf[5],pdf[6],pdf[7],pdf[8],pdf[9],pdf[10],pdf[11]) 
        end 
      end
     write(f,"---\n")
    end
  end
end

main()
