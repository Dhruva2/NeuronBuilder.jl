module NeuronBuilder

using ModelingToolkit 

const t = Num(Sym{Real}(:t))
const D = Differential(t)

include("channel_types.jl")
include("liu_channels.jl")
include("prinz_channels.jl")



include("synapses.jl")
include("build_neuron.jl")
include("compartments.jl")



export Synapse, IonChannel, Conductance, Soma, Compartment
export build_channel, build_synapse, build_neuron, build_group
export add_connection
export Liu_NaV, Liu_CaS, Liu_CaT, Liu_H, Liu_Ka, Liu_KCa, Liu_Leak, Liu_Kdr
export Prinz_NaV, Prinz_CaS, Prinz_CaT, Prinz_H, Prinz_Ka, Prinz_KCa, Prinz_Leak, Prinz_Kdr

# include("assets/Liu/Liu.jl")
# include("assets/Prinz/Prinz.jl")


end
