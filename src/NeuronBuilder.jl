module NeuronBuilder

using ModelingToolkit
using Unitful

const t = Num(Sym{Real}(:t))
const D = Differential(t)
export t, D

include("units.jl")
export Cm, Area, ICS, Params, Reversals

include("typetree.jl")
export Species, SpeciesProperty, Ion, SpeciesDynamics, Voltage, Sodium, Potassium, Calcium, Proton, Leak, Reversal, Current, Conductance, mRNA
export Compartment, Component, FlowChannel, Neuron, Synapse, EmptyConnection, ComponentSystem
export PlasticityRule, PlasticisedChannel
export Geometry, NoGeometry, capacitance

include("helper_functions.jl")
export get_name, shorthand_name, sensed, actuated, get_sensor, get_actuator, sensed_ions, voltage, sensedvars, vardivide
export reversals, currents, conductances, instantiate_variables, instantiate_parameters, instantiate_hooks

include("assets/UnitLiu/UnitLiu.jl")
using .UnitLiu
export UnitLiu

include("assets/Liu/Liu.jl")
using .Liu
export Liu

include("assets/Prinz/Prinz.jl")
using .Prinz
export Prinz

include("neurons.jl")
export BasicNeuron, BasicVoltageDynamics, EmptyNeuron

include("synapses.jl")
export Chol, Glut, directed_synapse

include("build_network.jl")
export build_network, build_group, add_connection, add_all_connections

include("plasticity.jl")
export OLeary_reg, OLearyCalcRegulation


end
