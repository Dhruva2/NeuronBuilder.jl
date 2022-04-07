1. Consistent way of turning parameters into states with dynamics. Regulated is one example. More generally want it to work to add synapses with dynamics.

2. Easy way of constructing networks of neurons. Maybe make an intermediate specification type. Will make life easier.

build_network(node_func, edge_func, args...; kwargs...)

node_func(i) = Neuron(i)
edge_func(i,j) = Synapse(i,j)


3. Find a nice way of adding callbacks. need to access e.g. the index of the voltage of the subsystems of different neurons in the solution array

VoltageIndex(ODESystem, i)  = index of voltage of ith neuron


4. ask andrea to get rid of the 14.96 and the prinz_conversion etc in the package. Just make an overall constant OUTSIDE of the package that multiplies the conductances. Will need different constants for calcium and voltage?