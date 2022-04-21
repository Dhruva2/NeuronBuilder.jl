### Basic components and subtypes ###

abstract type Species end
abstract type SpeciesProperty{S<:Species} end
abstract type Ion <: Species end
abstract type SpeciesDynamics{F} end

"""
Not including voltage in tagging ion channels. Flows are all charged particles. So they will always include voltage.
"""

struct Voltage <: Species end
struct Sodium <: Ion end
struct Potassium <: Ion end
struct Calcium <: Ion end
struct Proton <: Ion end
struct Chloride <: Ion end
struct MixedIon <: Ion end


abstract type Reversal{I<:Ion} <: SpeciesProperty{I} end
abstract type Current{I<:Ion} <: SpeciesProperty{I} end
abstract type Conductance{I<:Ion} <: SpeciesProperty{I} end
abstract type mRNA{I<:Ion} <: SpeciesProperty{I} end




abstract type Component end
abstract type Compartment <: Component end
abstract type Channel <: Component end

abstract type FlowChannel{Sensors<:Tuple,Actuators<:Tuple} <: Channel end
abstract type Synapse{Sensors_Pre,Sensors_Post,Actuators} <: Channel end

# channels always sense the reversal of the currents they actuate
FlowChannel(T) = FlowChannel{Tuple{Reversal{T}},Tuple{T}}
FlowChannel(S, A) = FlowChannel{Tuple{S,Reversal{A}},Tuple{A}}


# sensed from pre neuron, post neuron and actuated (onto post neuron)
Synapse(Spre, Spost, A) = Synapse{Tuple{Spre},Tuple{Spost},Tuple{A}}
struct EmptyConnection <: Synapse(Nothing, Nothing, Nothing) end


abstract type Neuron <: Compartment end
dynamics(n::Neuron) = n.dynamics
has_dynamics(n::Neuron, species) = haskey(dynamics(n), species)

# defined for reversal and conductance since they belong to a neuron. not conductance (belongs to Link)
has_dynamics(n::Neuron, ::Type{Current{I}}) where {I} = true
has_dynamics(n::Neuron, ::Type{Reversal{I}}) where {I} = haskey(dynamics(n), I)




abstract type Geometry end

struct NoGeometry{C} <: Geometry
    capacitance::C
end
capacitance(g::NoGeometry) = g.capacitance

abstract type PlasticityRule{S} end

struct PlasticisedChannel{S,S2,A} <: FlowChannel{Tuple{S,S2},A}
    channel::FlowChannel{S,A}
    modification::PlasticityRule{S2}
end

######### Geometries ###########
# is_geometric(::Compartment{T}) = T
# is_geometric(::Component) = false
"""
Every geometry needs methods for:
    capacitance()
    maybe reversals as dictionary?
    initial_conditions(f::Species)
"""


###### Dynamics #######
"""
Each flow dynamics struct needs a functor: 
    (s::SomeDynamics)(n::Neuron, variable, currents)

In the very near future need to add inputs and synaptic currents to this way of building
"""