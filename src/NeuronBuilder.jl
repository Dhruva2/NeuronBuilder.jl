module NeuronBuilder

using ModelingToolkit 

const t = Num(Sym{Real}(:t))
const D = Differential(t)

include("channel_types.jl")
include("synapses.jl")
include("build_neuron.jl")
include("compartments.jl")



export Synapse, IonChannel, Conductance, Soma, Compartment
export build_channel, build_synapse, build_neuron, build_group
export add_connection

include("assets/Liu/Liu.jl")
include("assets/Prinz/Prinz.jl")


end
