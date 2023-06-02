module NeuronBuilder

using ModelingToolkit

const t = Num(Sym{Real}(:t))
const D = Differential(t)
export t, D

include("Species.jl")
export Species, SpeciesProperty, Ion, SpeciesDynamics, Voltage, Sodium, Potassium, Calcium, Proton, Leak, Reversal, Current, Conductance, mRNA
export shorthand_name

include("Components.jl")
export Component, Compartment, FlowChannel, Neuron, Synapse, EmptyConnection, PlasticityRule, PlasticisedChannel, Geometry
export get_name

include("FlowChannels.jl")
export sensed, actuated, currents, sensedvars, mrna, reversals, conductances, sensed_ions, instantiate_variables, instantiate_parameters, instantiate_hooks, get_sensor, get_actuator, voltage

include("PlasticityRules.jl")
export OLeary_reg, OLearyCalcRegulation

include("SpeciesDynamics.jl")
export BasicVoltageDynamics
export get_parameters, get_states, default_params, default_states

include("Geometries.jl")
export NoGeometry
export capacitance

include("helper_functions.jl")
export vardivide


include("assets/Liu/Liu.jl")
using .Liu
export Liu

include("assets/Prinz/Prinz.jl")
using .Prinz
export Prinz


include("assets/HodgkinHuxley/HodgkinHuxley.jl")
using .HodgkinHuxley
export HodgkinHuxley

function channel_dynamics() end

include("neurons.jl")
export BasicNeuron, BasicVoltageDynamics, ResetDynamics, EmptyNeuron

include("synapses.jl")
export Chol, Glut, directed_synapse

include("Networks.jl")
export build_network, build_group, add_connection, add_all_connections

end
