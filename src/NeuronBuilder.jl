module NeuronBuilder

using ModelingToolkit

const t = Num(Sym{Real}(:t))
const D = Differential(t)
export t, D

include("typetree.jl")
export Species, Ion, SpeciesDynamics, Voltage, Sodium, Potassium, Calcium, Proton, Leak, Reversal
export Compartment, Component, FlowChannel, Neuron, Synapse, EmptyConnection, ComponentSystem
export PlasticityRule, PlasticisedChannel
export Geometry, NoGeometry, capacitance

include("helper_functions.jl")
export get_name, shorthand_name, sensed, actuated, get_sensor, get_actuator
export sensedvars, reversals, currents, conductance, instantiate_variables, instantiate_parameters

include("assets/Liu/Liu.jl")
using .Liu
export Liu

include("assets/Prinz/Prinz.jl")
using .Prinz
export Prinz

include("neurons.jl")
export BasicNeuron, BasicVoltageDynamics

include("synapses.jl")
export Chol, Glut

include("build_network.jl")
export build_network

include("plasticity.jl")
export OLearyCalcRegulation


end
