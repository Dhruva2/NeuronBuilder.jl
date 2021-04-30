module NeuronBuilder

using ModelingToolkit 

const t = Num(Sym{Real}(:t))
const D = Differential(t)

include("channels.jl")
include("synapses.jl")
include("build_neuron.jl")
include("compartments.jl")

export Synapse, IonChannel, Conductance, Soma, Compartment
export NaV, CaS, CaT, H, Ka, KCa, Kdr, Leak
export Chol, Glut
export build_channel, build_synapse, build_neuron, build_group
export add_connection

end
