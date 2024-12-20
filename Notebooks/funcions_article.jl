#### Funcions necessaries per a l'article exclúsivament ############

import LightGraphs as lg

"""
    create_qft_circuit_bis(n::Integer)

Function to create a qft circuit with n qubits
"""
function create_qft_circuit_bis(n::Integer)
    circ = QXZoo.Circuit.Circ(n)
    QXZoo.QFT.apply_qft!(circ, collect(1:n))
    QXZoo.Circuit.add_gatecall!(circ, QXZoo.DefaultGates.c_x(1, 2))
    QXZoo.Circuit.add_gatecall!(circ, QXZoo.DefaultGates.c_x(2, 3))
    circ
end


# **************************************************************************************** #
#                          Functions for creating various graphs
# **************************************************************************************** #

"""
    line_graph_tris(G::LabeledGraph)

Return a LabeledGraph representing the line graph of the 'G'. 

The label for each each vertex of the line graph is created by concatenating the labels of 
the corresponding vertices in the LabeledGraph 'G' but ordered
"""
function line_graph_tris(G::LabeledGraph)
    vertex_labels::Array{Symbol, 1}=[]
    
    G_edges = collect(QXGraphDecompositions.edges(G))
    #vertex_labels = [QXGraphDecompositions.combine_labels(G.labels[e.src], G.labels[e.dst]) for e in G_edges]
   
    for e in G_edges
       if e.src < e.dst
                etiqueta=Symbol(G.labels[e.src], :_, G.labels[e.dst])
                push!(vertex_labels,etiqueta)
                
            else
                etiqueta=Symbol(G.labels[e.dst], :_, G.labels[e.src])
                push!(vertex_labels,etiqueta)
            end
    end
    
    line_graph_tris(G.graph; vertex_labels=vertex_labels)
end

"""
    line_graph_tris(G::AbstractGraph;
               vertex_labels::Array{Symbol, 1}=Symbol[])

Return a LabeledGraph representing the line graph of the 'G'. 

The label for each each vertex of the line graph is created by 
concatenating the labels of the corresponding vertices in 'G'.

The symbols in the array `vertex_labels` are used as labels for the vertices of the returned
line graph. If `vertex_labels` is empty then labels are created by combining the indices of 
the corresponding vertices in 'G'.
"""
function line_graph_tris(G::lg.AbstractGraph; 
                    vertex_labels::Array{Symbol, 1}=Symbol[])
    # Create a labeled graph LG whose vertices corresponding to the edges of G.
    
    G_edges = collect(lg.edges(G))
    if isempty(vertex_labels)
        #vertex_labels = [QXGraphDecompositions.combine_labels(e.src, e.dst) for e in G_edges]
        #vertex_labels=[]
        for e in G_edges
            if e.src < e.dst
                etiqueta=Symbol(G.labels[e.src], :_, G.labels[e.dst])
                push!(vertex_labels,etiqueta)
                
            else
                etiqueta=Symbol(G.labels[e.dst], :_, G.labels[e.src])
                push!(vertex_labels,etiqueta)
            end
        end
    end
    LG = LabeledGraph(lg.SimpleGraph(length(G_edges)), vertex_labels)

    # Connect any two vertices of LG whose corresponding edges in G share a vertex in G.
    for i in 1:length(G_edges)-1
        for j in i+1:length(G_edges)
            u = G_edges[i]
            v = G_edges[j]
            if (u.src == v.src) || (u.src == v.dst) || (u.dst == v.src) || (u.dst == v.dst)
                QXGraphDecompositions.add_edge!(LG, i, j)
            end
        end
    end
    LG
end



"""
   GN_pla(G::AbstractGraph,;
               tnc_tn::TensorNetwork)

Return a  contraction plan based on Girvan-Newman Algorith (GN). 


"""


 
function GN_pla(g,tnc_tn)
        
     h_ig=lg2ig(g) 

     comunitat=ordenacio_girvan_igraph(h_ig)
     
     ordenacio=Vector{Symbol}[]
     for i in comunitat
        parella=[]
    
        tensor_1= Symbol("t"*(string(i[1])))
        tensor_2= Symbol("t"*(repr(i[2]))) 
        push!(parella,tensor_1,tensor_2)
        push!(ordenacio,parella)
     end
      
        pla=order_to_contraction_plan(ordenacio,tnc_tn)
  end


function GN_pla(g,vmap,tnc_tn)
        
     h_ig=lg2ig(g) 

     comunitat=ordenacio_girvan_igraph(h_ig)
     
     ordenacio=Vector{Symbol}[]
     for i in comunitat
        parella=[]
    
        tensor_1= Symbol("t"*(string(vmap[i[1]])))
        tensor_2= Symbol("t"*(repr(vmap[i[2]]))) 
        push!(parella,tensor_1,tensor_2)
        push!(ordenacio,parella)
     end
      
        pla=order_to_contraction_plan(ordenacio,tnc_tn)
  end




# **************************************************************************************** #
#       Converting edge elimination orders & contraction trees into contraction plans
# **************************************************************************************** #

"""
    order_to_contraction_plan(elimination_order::Array{<:Array{Symbol, 1}, 1},
                              tn::Union{TensorNetwork, OrderedDict{Symbol, Array{Index, 1}}}
                              )::Array{NTuple{3, Symbol}, 1}
Convert the given edge elimination order into a contraction plan for `tn`.
`tn` can be a TensorNetwork or an OrderedDict describing a set of tensors, in which case
the keys of `tn` are assumed to be tensor ids/names and the values are arrays containing the
indices of the corresponding tensor.
"""
function order_to_contraction_plan(elimination_order::Array{Array{Symbol, 1}, 1},
                                   tn::Union{TensorNetwork, OrderedDict{Symbol, Array{Index, 1}}}
                                   )::Array{NTuple{3, Symbol}, 1}
    # An array to hold the contraction plan.
    plan = Array{NTuple{3, Symbol}, 1}()

    # A dictionary to keep track of which tensors are replaced by which intermediate tensors
    # at different stages of the contraction process. Initially, before any pairwise
    # contractions, none of the tensors are replaced by intermediates, so all tensor ids are
    # mapped to themselves.
    intermediate_tensor_map = Dict{Symbol, Symbol}(keys(tn) .=> keys(tn))

    # Convert each edge in the elimination order to a set of pairwise contractions and
    # append it to plan.
    for (i, edge) in enumerate(elimination_order)
        if length(edge) == 2
            # For edges consisting of just two tensors, find the intermediates they belong
            # to and add their pairwise contraction to the plan.
            A_id = _get_intermediate_tensor(intermediate_tensor_map, edge[1])
            B_id = _get_intermediate_tensor(intermediate_tensor_map, edge[2])
            if !(A_id == B_id)
                # id for the new intermediate created by contracting A_id and B_id.
                I_id = Symbol("I$i")
                append!(plan, [(A_id, B_id, I_id)])

                # Update intermediate_tensor_map with new intermediate.
                intermediate_tensor_map[A_id] = I_id
                intermediate_tensor_map[B_id] = I_id
                intermediate_tensor_map[I_id] = I_id
            end

        elseif length(edge) > 2
            # For edges with more than 2 tensors, collect all of the tensors and
            # intermediate tensors that belong to the edge and find a contraction plan for
            # them. Append the contraction this plan to plan.
            tensors_to_contract = OrderedDict{Symbol, Array{Index, 1}}()
            for t_id in edge
                I_id = _get_intermediate_tensor(intermediate_tensor_map, t_id)
                inds = typeof(tn) <: TensorNetwork ? QXTns.inds(tn[t_id]) : tn[t_id]
                tensors_to_contract[I_id] = symdiff(get(tensors_to_contract, I_id, []), inds)
            end
            if length(tensors_to_contract) > 1
                # Check if netcon can be used on the given set of tensors. If not, use
                # a fallback method to find a contraction plan.
                tensor_sizes = prod.([QXTns.dim.(inds) for inds in values(tensors_to_contract)])
                if length(tensors_to_contract) < 37 && any(tensor_sizes .< 2^62)
                    local_contraction_plan = netcon(tn, tensors_to_contract)
                else
                    local_contraction_plan = min_fill_contraction_plan(tensors_to_contract)
                end
                append!(plan, local_contraction_plan)

                # Update intermediate_tensor_map with new intermediates.
                for (A_id, B_id, I_id) in local_contraction_plan
                    intermediate_tensor_map[A_id] = I_id
                    intermediate_tensor_map[B_id] = I_id
                    intermediate_tensor_map[I_id] = I_id
                end
            end
        end
    end
    plan
end




"""
    _get_intermediate_tensor(intermediate_tensor_map::Dict{Symbol, Symbol}, t_id::Symbol)::Symbol
Return the id of the intermediate tensor of which the tensor `t_id` is a factor.
During a tensor network contraction, tensors and intermediate tensors are contracted
together and replaced by the result. The dictionary `intermediate_tensor_map` is
assumed to describe which intermediate tensor has replaced a given tensor during a
contraction of a tensor network.
"""
function _get_intermediate_tensor(intermediate_tensor_map::Dict{Symbol, Symbol}, t_id::Symbol)::Symbol
    while true
        I_id = intermediate_tensor_map[t_id]
        # If t_id = I_id then t_id has not been replaced by an intermediate tensor yet.
        (t_id == I_id) && return I_id
        t_id = I_id
    end
end



"""
    lg2ig(lg_g)(G::AbstractGraph)
 Convert a Julia graph into a Python graph. See 

https://python.igraph.org/en/latest/index.html
"""

ig = pyimport("igraph");

function lg2ig(lg_g)
           i_g = ig.Graph(LightGraphs.nv(lg_g))
           arestes=[]
           for e in LightGraphs.edges(lg_g)
               a=(LightGraphs.src(e)-1,LightGraphs.dst(e)-1)
               #println(a)
               push!(arestes,a)
           end
           i_g.add_edges(arestes)
           return i_g
       end
########################################



"""
ordenacio_girvan_igraph(h)
GN algorith in Python graph. See 

https://python.igraph.org/en/latest/index.html
"""

#Fem una versió del algoritme de GN per a igraph
function ordenacio_girvan_igraph(h)
    eixida=[]
    
    for i in 1:length(h.get_edgelist())

        edge_betweenness = h.edge_betweenness()
        aresta=findmax(edge_betweenness)[2]
        e=h.es[aresta]
        edge = e.source+1, e.target+1 
        push!(eixida,edge)
        id=h.get_eid(e.source, e.target) 
        h.delete_edges(id)

    end
    eixida_inv=reverse(eixida)
    return eixida_inv
end


"""
labelg_to_communitats_between(labeled__light_graf,clusters)

This function gives us betweenness communities and its modularity
See 

https://python.igraph.org/en/latest/index.html

"""


function labelg_to_communitats_between(labeled__light_graf,clusters)
    
        h_Lg_ig=lg2ig(labeled__light_graf.graph);
        if clusters==0
            communities_Lg =h_Lg_ig.community_edge_betweenness( ) 
        else
            communities_Lg =h_Lg_ig.community_edge_betweenness(clusters=clusters) 
        end
    
        communities_Lg = communities_Lg.as_clustering();
        comunitats_betwenness=py"""list"""(communities_Lg);
        comunitats_julia=[];
        for i in comunitats_betwenness
            planet=Int[]
            for j in 1:length(i)
            push!(planet,i[j]+1)
            end
        push!(comunitats_julia,planet)
        end
        modularitat=h_Lg_ig.modularity(communities_Lg) # aquest pas no es necessari dins la funció si he de promitjar
    return comunitats_julia,comunitats_betwenness,modularitat
end

"""
pla_contraccio_multiple_G_N(nova_llista_comunitats_julia,tnc, g)


We obtain a contraction plan for each community

"""



#########################################################################
#Obtenim un pla de contracció per a cada comunitat

import QXGraphDecompositions

function pla_contraccio_multiple_G_N(nova_llista_comunitats_julia,tnc, g)

plans=[]
for i in 1:length(nova_llista_comunitats_julia)
   subgraf=nova_llista_comunitats_julia[i]
   sg, vmap = LightGraphs.induced_subgraph(g, subgraf)
   vmap_s =Symbol.(vmap)
         
  pla_comunitat=GN_pla(sg,vmap,tnc.tn)
  push!(plans,pla_comunitat)
end
#println(plans)
    #############################################################
    
    
    # passem les comunitats numériques a tensors 
comunitats_x=[]
for j in 1:length(nova_llista_comunitats_julia)
   comunitat=[ Symbol("t$i")   for i in nova_llista_comunitats_julia[j]]
   push!(comunitats_x,comunitat) 
end
#println(comunitats_x)

# En funció del nombre de comunitats, creem un vector de xarxes tensorials
tns=[]
 for j in 1:length(nova_llista_comunitats_julia)
    tensors_c=[]
    for i in comunitats_x[j]
       push!(tensors_c,tnc.tn.tensor_map[i])
    end
    tensors_c=Vector{ QXTensor}(tensors_c)
  #println(tensors_c3)
    tn=TensorNetwork(tensors_c,comunitats_x[j])
    push!(tns,tn)
end

   return tns,plans
end





using ITensors
using LinearAlgebra
using NDTensors

"""
    TensorNetwork(array::Vector{<: QXTensor})

Outer constructor to create a tensor network object (subnetwork) from an array of ITensor objects and an array of symbols as names
"""
function TensorNetwork(array::Vector{<: QXTensor},comunitat::Vector{Symbol})
    tensor_map = OrderedDict{Symbol, QXTensor}()
    bond_map = OrderedDict{Index, Vector{Symbol}}()
    next_id = 1
    for t in array
        tensor_id = Symbol(comunitat[next_id])
        next_id += 1
        tensor_map[tensor_id] = t
        for bond in inds(t)
            if haskey(bond_map, bond)
                push!(bond_map[bond], tensor_id)
            else
                bond_map[bond] = [tensor_id]
            end
        end
    end
    TensorNetwork(tensor_map, bond_map, next_id)
end


"""
    primera_contraccio_paral(tns,plans)

These functions have as input the vectors of TNS circuits and the planes and give us the unique network  for the final contraction
"""

# aquestes funcions tenen com a entrada els vectors de circuits tns i els plans i ens donen la xarxa única c
function primera_contraccio_paral(tns,plans)
    contraccio_paral(tns, plans)

    c=TensorNetwork()
    for i in 1: length(tns)
   
          c=Base.merge(c, tns[i])
    end

    return c 
   
end


using DataStructures
using Distributed

function contraccio_paral(tns,plans)
   
    
   @sync   Base.Threads.@threads for i in 1:length(tns)
       contract_tn!(tns[i], plans[i])
    end
     tns
end


"""
    min_fill_contraction_plan_tw(tn::TensorNetwork;
                                   hypergraph::Bool=false)

this function returns a plan using min_fill algorith and its treewidth
"""



import QXGraphDecompositions as qxg
import LightGraphs

function min_fill_contraction_plan_tw(tn::TensorNetwork;
                                   hypergraph::Bool=false)
    # Convert tn to a line graph and pass it to flow cutter to find an tree decomposition.
    lg, symbol_map = convert_to_line_graph(tn, use_hyperedges=hypergraph)

    # Use flow cutter to try find a tree decomposition of the line graph.
    tw, order = qxg.min_fill(lg)

    # Convert the elimination order to an array of Index structs in tn.
    order = [symbol_map[lg_vertex] for lg_vertex in order]

    # Convert the elimination order into a contraction plan.
    pla = order_to_contraction_plan(order, tn)
    return tw,pla 
end

min_fill_contraction_plan_tw(tnc::TensorNetworkCircuit; kwargs...) = min_fill_contraction_plan_tw(tnc.tn; kwargs...)


"""
    pla_paral_p(c,pla)

This function returns an ordered plane with the parallelizable elements - 
It presents as inputs the final single c network after the first parallel contraction and the sequential plane which can be minfill, FlowCutter or other.
"""


##########################################################################
# Aquesta funció retorna un pla ordenat amb l elements paral·lelitzables 
# Presenta com a entrades la xarxa c única i el pla seqüencial que pot ser minfill,Fc o cap altre.
#################################################################################
function pla_paral_p(c,pla)
    b_p=[];b_q=[];

   for i in pla
    #println(i)
     if i[1] in keys(c) && i[2] in keys(c)
        push!(b_p,i)
     else
       push!(b_q,i) 
     end
    
    end
         l= length(b_p)
         pla_p = append!(b_p,b_q)
    return pla_p,l
end

"""
    contrau_p(c_gn,pla_mf_p,p)

This function combined with the previous one serves to contract with a certain parallelism the third phase of the algorithm or phase of the final contraction. Returns the probability amplitude that a given input give a given output
"""

 #################################################################################


function contrau_p(c_gn,pla_mf_p,p)
                 @sync   Base.Threads.@threads for i in 1:p
                      
                      contract_pair!(c_gn, pla_mf_p[i][1],pla_mf_p[i][2],pla_mf_p[i][3])
                  end
                for j in p+1:length(pla_mf_p)
                        contract_pair!(c_gn, pla_mf_p[j][1],pla_mf_p[j][2],pla_mf_p[j][3])
                  end
    return c_gn.tensor_map[pla_mf_p[end][3]].storage
        end

#############################################################

"""
    Calcul_GN_Sequencial(cct,timings)

Sequential contraction of the tensor network using one of the results of the Girvan–Newman algorithm as a contraction plan
"""

using TimerOutputs

function Calcul_GN_Sequencial(cct,timings)

   
  if timings
        reset_timer!()
    end

  @timeit " 1T.Obtaining a line graph" begin


     tnc= convert_to_tnc(cct)

    light_graf = convert_to_graph(tnc)
    
    labeled__light_graf=LabeledGraph(light_graf);

    labeled_line_light_graf= line_graph_tris(labeled__light_graf); # obtenim el line graf del graf original

  end
 
  @timeit " 2T.Getting GN plan" begin

      pla_gn=GN_pla(light_graf,tnc.tn)

   end
  
  @timeit " 3T.Final contraction " begin
     s=contract_tn!(tnc.tn, pla_gn)

  end

 if timings
            print_timer()
        end

return s

end


"""
   ComParCPU(circ,entrada,eixida,n_com;timings=true,decompose=true)

This function would be our algorithm in 3 phases described in the paper. 

Inputs: 
circ is the circuit we want to evaluate
entrada means input qubits
eixida means output qubits
n_com is the number of communities you want to explore
timings: true means graphical results
decompose: if true, cicuit gates are decomposed.

Output: s as a solution that represents an amplitude of probability
"""


using TimerOutputs
#using CUDA

# Aquesta funció seria el nostre algorisme en 3 fasses on podem modular: algorisme primera contraccio entre gn i minfill i
# la segona contracció en sèrie (només minfill, però podriem ampliar-ho) o en paral·lel

function ComParCPU(circ,entrada,eixida,n_com;timings=true,decompose=true)

    if timings
        reset_timer!()
    end



          @timeit "1T.Obtaining Communities" begin
               num_qubits = circ.num_qubits              
               tnc= convert_to_tnc(circ;no_input=false,no_output=false,input=entrada,output=eixida,decompose=true )                 
              
               light_graf  = convert_to_graph(tnc)



         
               labeled__light_graf=LabeledGraph(light_graf);
               labeled_line_light_graf= line_graph_tris(labeled__light_graf); 
   
                h_Lg_ig=lg2ig(labeled__light_graf.graph);
                h_Lg_ig.summary(verbosity=1);



                 comunitats_julia,comunitats_betwenness,modularitat=labelg_to_communitats_between(labeled__light_graf,n_com);
  
           end 
    
    
    
            @timeit "2T.Parallel contraction of communities" begin
  
  
                    tns,plans= pla_contraccio_multiple_G_N(comunitats_julia,tnc,light_graf)
                    c =primera_contraccio_paral(tns,plans);
           
            end 
      
    

 


                @timeit "3T.Final contraction" begin
              

              tw,pla= min_fill_contraction_plan_tw(c)
               
                s=contract_tn!(c, pla)
                
             end

    

    if timings
            print_timer()
        end

return s
end

"""
   ComParCPU_para(circ,entrada,eixida,n_com;timings=true,decompose=true)

This function would be our algorithm in 3 phases described in the paper with a modification in the final contraction where we can
 parallelize the last contraction. 

Inputs: 
circ is the circuit we want to evaluate
entrada means input qubits
eixida means output qubits
n_com is the number of communities you want to explore
timings: true means graphical results
decompose: if true, cicuit gates are decomposed.

Output: s as a solution that represents an amplitude of probability
"""



using TimerOutputs
#using CUDA

# Aquesta funció seria el nostre algorisme en 3 fasses on podem modular: algorisme primera contraccio entre gn i minfill i
# la segona contracció en  paral·lel, com aquest cas o en sèrie com abans.

function ComParCPU_para(circ,entrada,eixida,n_com;timings=true,decompose=true)

    if timings
        reset_timer!()
    end



          @timeit "1T.Obtaining Communities" begin
               num_qubits = circ.num_qubits              
               tnc= convert_to_tnc(circ;no_input=false,no_output=false,input=entrada,output=eixida,decompose=true )                 
              
               light_graf  = convert_to_graph(tnc)



         
               labeled__light_graf=LabeledGraph(light_graf);
               labeled_line_light_graf= line_graph_tris(labeled__light_graf); 
   
                h_Lg_ig=lg2ig(labeled__light_graf.graph);
                h_Lg_ig.summary(verbosity=1);



                 comunitats_julia,comunitats_betwenness,modularitat=labelg_to_communitats_between(labeled__light_graf,n_com);
  
           end 
    
    
    
            @timeit "2T.Parallel contraction of communities" begin
  
  
                    tns,plans= pla_contraccio_multiple_G_N(comunitats_julia,tnc,light_graf)
                    c =primera_contraccio_paral(tns,plans);
           
            end 
      
    

 


               
                    @timeit "3T.Final contraction in parallel" begin
                
                      
                          tw,pla= min_fill_contraction_plan_tw(c)
                          pla_mf_p,p  = pla_paral_p(c,pla)
                          s= contrau_p(c,pla_mf_p,p)
                    end

    if timings
            print_timer()
        end

return s
end
