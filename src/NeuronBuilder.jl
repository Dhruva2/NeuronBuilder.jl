module NeuronBuilder

# Write your package code here.
using ModelingToolkit, OrdinaryDiffEq, Unitful

include("channels.jl")
include("synapses.jl")
include("build_neuron.jl")

export Synapse, IonChannel, Conductance
export NaV, CaS, CaT, H, Ka, KCa, Kdr, Leak
export Chol, Glut
export build_channel, build_synapse, build_neuron, build_group
export add_connection


end
