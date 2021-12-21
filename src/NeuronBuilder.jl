module NeuronBuilder

using ModelingToolkit

const t = Num(Sym{Real}(:t))
const D = Differential(t)
export t, D

include("channel_types.jl")
export Synapse, IonChannel, Soma, Compartment, ComponentSystem

function channel_dynamics() end
function ionic_current() end
function calcium_current() end
function external_params() end

include("assets/Liu/Liu.jl")
using .Liu
export Liu

include("assets/Prinz/Prinz.jl")
using .Prinz
export Prinz

include("synapses.jl")
include("helper_functions.jl")
include("build_neuron.jl")

export channel_dynamics

export build_channel, build_neuron, build_group, add_connection

end
