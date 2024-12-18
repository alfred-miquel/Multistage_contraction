

using QXTools;
using Graphs;
using QXTns;
using QXZoo;
using PyCall;


include("introduccio.jl");

include("contraccions.jl");

include("contraccions_bis.jl");


#########################################
#per als minfill de la pace 2017
include("algorismes.jl")
include("minfill_pace_2017.jl")

###########################################


@info("Quants qubits vols que proven en els algorismes (n)? \n\n")

# Calling rdeadline() function
N = readline()
n = parse(Int64, N)
# Creacio del circuit

@info("Creació del circuit cct")

#n=20;

println("El circuit QFT és de $n bits")
cct=create_qft_circuit_bis(n);


#------------------------------------------------------------



### Fem una funció per a posar timers############

 using TimerOutputs

function Calcul_GN_Sequencial(cct,timings)

   
  if timings
        reset_timer!()
    end

  @timeit " 1T.Obtenció de line graf" begin


     tnc= convert_to_tnc(cct)

    light_graf = convert_to_graph(tnc)
    
    labeled__light_graf=LabeledGraph(light_graf);

    labeled_line_light_graf= line_graph_tris(labeled__light_graf); # obtenim el line graf del graf original

  end
 
  @timeit " 2T.Obtenció del pla GN" begin

      pla_gn=GN_pla(g,tnc.tn)

   end
  
  @timeit " 3T.Contracció sequencial usant BLAS" begin
     s=contract_tn!(tnc.tn, pla_gn)

  end

 if timings
            print_timer()
        end

return s

end

