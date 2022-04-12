module NeuronBuilder

using ModelingToolkit
using Unitful

const t = Num(Sym{Real}(:t))
const D = Differential(t)
export t, D

include("units.jl")
export Cm, Area, ICS, Params, Reversals

include("channel_types.jl")
export Synapse, IonChannel, RegIonChannel, RegIon, Soma, Compartment, Component, ComponentSystem, Liu_conversion, Prinz_conversion, syn_conv_factor, EmptyConnection

include("specification_types.jl")
export Neuron
include("build_network.jl")
export build_network

include("helper_functions.jl")
export get_name, get_g, Area, Cm, L2NB, P2NB
function channel_dynamics() end
function ionic_current() end
function calcium_current() end
function external_params() end

include("assets/UnitLiu/UnitLiu.jl")
using .UnitLiu
export UnitLiu

include("assets/Liu/Liu.jl")
using .Liu
export Liu

include("assets/Prinz/Prinz.jl")
using .Prinz
export Prinz

include("synapses.jl")
include("regulation.jl")
export Regulated

include("build_neuron.jl")

export channel_dynamics

export build_channel, build_neuron, build_group, add_connection
export Chol, Glut

end
