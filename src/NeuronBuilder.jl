module NeuronBuilder

# Write your package code here.
using ModelingToolkit, OrdinaryDiffEq

include("initialisations.jl")
include("channels.jl")
include("build_neuron.jl")

export VoltageChannel, CalciumChannel
export NaV, CaS, CaT, H, Ka, KCa, Kdr, Leak
export Synapse, Chol, Glut
export build_channel, build_neuron
export voltage_current, calcium_current

end
