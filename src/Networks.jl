"""
What to do:
"""

"""
nodefn(i) = neuron at node i on graph
edgejn(i,j) = synapse connecting nodes i and j

- remember that neuron(j) needs to be instantiated with all synapses ∑_i node(i,j) in its preceding field
- need a hook() function that connects all variables of type Previous

for i, j
    build list of connecting equations using hook()
    build list of subsystems (i.e. neurons )
end

return ODESystem(connecting+equations; systems = all_subsystems)
"""
function build_network(nodefn, edgefn)
    
end

"""
1. synapse_ODESystem = getproperty(post_neuron, :synapse_name::Symbol)
2. for each variable of type Previous{Q<:Quantity}
        pre_neuron.Q ~ post_neuron.synapse.Previous{Q}
end
"""
function hook(pre_neuron, post_neuron, synapse_name::Symbol)
    
end