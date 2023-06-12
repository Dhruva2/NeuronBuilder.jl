module NeuronBuilder

using ModelingToolkit

const t = Num(Sym{Real}(:t))
const D = Differential(t)
export t, D

include("Species.jl")
export TrackedQuantity, UntrackedQuantity

export Species, SpeciesProperty, Ion, SpeciesDynamics, Voltage, Sodium, Potassium, Calcium, Proton, Leak, Reversal, Current, Conductance, mRNA, shorthand_name, Capacitance

export Quantity, UntrackedQuantity, NoGeometry

export Previous, Post, OrderRelation


include("Components.jl")
export Compartment, Component, BasicChannel, DirectedChannel, Neuron, Geometry, PlasticityRule, EmptySynapse
export build_vars, build_vars_from_owner, tagged_internal_variables, untagged_internal_variables, sensed, actuated, required_defaults

export find_from



function sensed() end
function actuated() end
function is_dynamic() end
function has_dynamics() end #delete this one prob
function tagged_internal_variables() end
function untagged_internal_variables() end 

include("../assets/BasicComponents/BasicComponents.jl")
using .BasicComponents
export BasicComponents



include("../assets/Liu/Liu.jl")
using .Liu
export Liu


include("Neurons.jl")
export BasicNeuron
export geometry, dynamics, defaults, build_defaults
end
