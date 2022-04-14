### Basic components and subtypes ###
abstract type Species end
abstract type Ion <: Species end
abstract type PseudoIon <: Ion end
abstract type SpeciesDynamics{F} end

"""
Not including voltage in tagging ion channels. Any flow channel. Flows are all charged particles. So they will always include voltage. Not going to consider flows of e.g. proteins
"""
struct Voltage <: Species end
struct Sodium <: Ion end
struct Potassium <: Ion end
struct Calcium <: Ion end
struct Proton <: Ion end
struct Leak <: PseudoIon end
abstract type Reversal{I<:Ion} end

ionic(x) = x <: Ion
export ionic

abstract type Component end
abstract type Compartment <: Component end
abstract type FlowChannel{Sensors,Actuators} <: Component end

abstract type Neuron <: Compartment end

abstract type Synapse <: Component end
struct EmptyConnection <: Synapse end

struct ComponentSystem{C<:Component,S<:AbstractTimeDependentSystem}
    c::C
    sys::S
end

### channels always sense the reversal of the currents they actuate
# generalize for more than actuator
FlowChannel(T) = FlowChannel{Reversal{T},T}
FlowChannel(S, A) = FlowChannel{Tuple{S,Reversal{A}},A}

abstract type Geometry end
struct NoGeometry{C} <: Geometry
    capacitance::C
end
capacitance(g::NoGeometry) = g.capacitance

abstract type PlasticityRule{S} end

function type_ps(t::Type)
    isempty(t.parameters) && return (t,)
    return t.parameters |> Tuple
end
function merge_types(t1::Type, t2::Type)
    ps = (type_ps(t1)..., type_ps(t2)...) |> unique
    return Tuple{ps...}
end

struct PlasticisedChannel{S,S2,A} <: FlowChannel{merge_types(Type{S}, Type{S2}),A}
    channel::FlowChannel{S,A}
    mutation::PlasticityRule{S2}
end


######### Geometries ##############
# is_geometric(::Compartment{T}) = T
# is_geometric(::Component) = false
"""
Every geometry needs methods for:
    capacitance()
    maybe reversals as dictionary?
    initial_conditions(f::Species)
"""


#### Dynamics #######
"""
Each flow dynamics struct needs a functor: 
    (s::SomeDynamics)(n::Neuron, variable, currents)

In the very near future need to add inputs and synaptic currents to this way of building
"""