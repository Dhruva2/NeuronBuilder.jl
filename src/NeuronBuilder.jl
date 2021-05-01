module NeuronBuilder

# Write your package code here.
using ModelingToolkit, OrdinaryDiffEq, Unitful

include("channel_types.jl")
include("synapses.jl")

    module Liu
    import NeuronBuilder
    include("assets/Liu/channels.jl")
    export NaV, CaS, CaT, Leak
    end



include("build_neuron.jl")
include("compartments.jl")

#include("/libraries/Liu/Liu.jl")

export Synapse, IonChannel, Conductance, Soma, Compartment
export build_channel, build_synapse, build_neuron, build_group
export add_connection




end
