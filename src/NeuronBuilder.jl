module NeuronBuilder

using ModelingToolkit

const t = Num(Sym{Real}(:t))
const D = Differential(t)
export t, D

include("channel_types.jl")
export Synapse, IonChannel, Conductance, Soma, Compartment

function channel_dynamics() end

include("assets/Liu/Liu.jl")
using .Liu
export Liu


include("assets/Prinz/Prinz.jl")
using .Prinz
export Prinz

export channel_dynamics

include("synapses.jl")
include("build_neuron.jl")
include("compartments.jl")


export build_channel, build_synapse, build_neuron, build_group
export add_connection
#export Liu_NaV, Liu_CaS, Liu_CaT, Liu_H, Liu_Ka, Liu_KCa, Liu_Leak, Liu_Kdr
#export Prinz_NaV, Prinz_CaS, Prinz_CaT, Prinz_H, Prinz_Ka, Prinz_KCa, Prinz_Leak, Prinz_Kdr

# 
# include("assets/Prinz/Prinz.jl")





end
