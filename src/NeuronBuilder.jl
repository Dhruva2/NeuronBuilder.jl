module NeuronBuilder

using ModelingToolkit

const t = Num(Sym{Real}(:t))
const D = Differential(t)
export t, D

include("typetree.jl")
export Synapse, IonChannel, RegIonChannel, RegIon, Soma, Compartment, Component, ComponentSystem, Liu_conversion, Prinz_conversion, syn_conv_factor, EmptyConnection

export Voltage, Ion, Sodium, Potassium, Calcium, Proton, Leak, Reversal, Species
export FlowChannel, FlowDynamics, BasicVoltageDynamics, LiuCalciumDynamics, LiuCaReversalDynamics
export ions

export Injection
export PlasticityRule

export ComponentInfo, IonicConductances, Flows, IonicCurrents, CalciumFlows, ExternalHooks, InstantiateVariables, Reversal


include("neurons.jl")
export Neuron, BasicNeuron
export Geometry, NoGeometry
export capacitance
include("build_network.jl")
export build_network

include("helper_functions.jl")
export get_name, get_g, Area, Cm, L2NB, P2NB
function channel_dynamics() end
function ionic_current() end
function ionic_conductance() end
function calcium_current() end
function external_params() end



include("assets/Liu/Liu.jl")
using .Liu
import .Liu
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
