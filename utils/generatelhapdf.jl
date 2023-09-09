#!/usr/bin/julia
using BAT, DensityInterface
using PartonDensity
using QCDNUM
using Random, Distributions, ValueShapes, ParallelProcessingTools
using StatsBase, LinearAlgebra
using Printf
using ArgParse
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
        "--replicas", "-r"
            help = "replicas"
            arg_type = Int
            default = 100
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

gg = QCDNUM.load_params(string("fitresults/", parsed_args["fitresults"], "2.h5"))
qcdnum_params =  gg["evolution_params"]
splint_params = QCDNUM.SPLINTParams();
quark_coeffs = QuarkCoefficients();
forward_model_init(qcdnum_params, splint_params)

allx = exp10.(range(-3.0,0.0,length=20))
allq = exp10.(range(1.0,2.0,length=20))
allalpha = copy(allq)
for (i, qq) in enumerate(allq)
  allalpha[i] =QCDNUM.asfunc(qq*qq)[1]
end

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
OrderQCD: ",string(qcdnum_params.order),"
ForcePositive: 1
FlavorScheme: fixed
NumFlavors: ",string(qcdnum_params.n_fixed_flav),"
XMin: ",string(qcdnum_params.grid_params.x_min[1]),"
XMax: 1
QMin: ",string(qcdnum_params.q0),"
QMax: ",string(qcdnum_params.grid_params.qq_bounds[2]),"
MZ: 91.1876
MUp: 0
MDown: 0
MStrange: 0
MCharm: 1.3
MBottom: 4.75
MTop: 172
AlphaS_MZ: ",string(qcdnum_params.α_S),"
AlphaS_OrderQCD: ",string(qcdnum_params.order),"
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

    samples_data = bat_read((string("fitresults/", parsed_args["fitresults"], ".h5"),"samples")).result;
    sub_samples_all = BAT.bat_sample(samples_data, BAT.OrderedResampling(nsamples=parsed_args["replicas"])).result;
    samples_mode2 = mode(sub_samples_all);
    sub_samples = BAT.bat_sample(samples_data, BAT.OrderedResampling(nsamples=Ns)).result;

    pdf_params_s = DirichletPDFParams(K_u=samples_mode2.K_u, K_d=samples_mode2.K_d, K_q=samples_mode2.K_q,
                                      λ_g1=samples_mode2.λ_g1, 
                                      λ_g2=samples_mode2.λ_g2,
                                      K_g=samples_mode2.K_g, 
                                      λ_q=samples_mode2.λ_q, 
                                      θ=Vector(samples_mode2.θ))
    allparams = [pdf_params_s]

    for s in eachindex(sub_samples)
      push!(allparams,DirichletPDFParams(K_u=sub_samples.v.K_u[s], K_d=sub_samples.v.K_d[s], K_q=sub_samples.v.K_q[s],
                                      λ_g1=sub_samples.v.λ_g1[s], 
                                      λ_g2=sub_samples.v.λ_g2[s],
                                      K_g=sub_samples.v.K_g[s], 
                                      λ_q=sub_samples.v.λ_q[s], 
                                      θ=Vector(sub_samples.v.θ[s]))
      )
    end

    NN=0
    #push!(allparams,samples_mode2)
    for s in  allparams
    
     my_func = get_input_pdf_func(s)
     input_pdf = QCDNUM.InputPDF(func=my_func, map=input_pdf_map)

    # Evolve PDFs
     iq0 = QCDNUM.iqfrmq(qcdnum_params.q0)
     ϵ = QCDNUM.evolfg(qcdnum_params.output_pdf_loc, input_pdf.cfunc, input_pdf.map, iq0)
     itype=1

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
          qq2=qq*qq
          pdf = QCDNUM.allfxq(itype, xx, qq2, 0, 1)
          @printf(f, "%8.8f %8.8f %8.8f %8.8f %8.8f %8.8f %8.8f %8.8f %8.8f %8.8f %8.8f\n",pdf[1],pdf[2],pdf[3],pdf[4],pdf[5],pdf[6],pdf[7],pdf[8],pdf[9],pdf[10],pdf[11]) 
        end 
      end
     write(f,"---\n")
    end
  end
end

main()
