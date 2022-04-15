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

shorthand_name(::Type{Voltage}) = :V
shorthand_name(::Type{Sodium}) = :Na
shorthand_name(::Type{Potassium}) = :K
shorthand_name(::Type{Calcium}) = :Ca
shorthand_name(::Type{Proton}) = :H
shorthand_name(::Type{Leak}) = :Leak
shorthand_name(::Type{Reversal{T}}) where {T} = Symbol(:E, shorthand_name(T))
shorthand_name(x::Type{Tuple{T,R}}) where {T,R} = shorthand_name.(x.types)


ionic(x) = x <: Ion
export ionic

abstract type Component end
abstract type Compartment <: Component end
abstract type FlowChannel{Sensors<:Tuple,Actuators<:Tuple} <: Component end

abstract type Neuron <: Compartment end

abstract type Synapse <: Component end
struct EmptyConnection <: Synapse end

struct ComponentSystem{C<:Component,S<:AbstractTimeDependentSystem}
    c::C
    sys::S
end

### channels always sense the reversal of the currents they actuate
# generalize for more than actuator
FlowChannel(T) = FlowChannel{Tuple{Reversal{T}},Tuple{T}}
FlowChannel(S, A) = FlowChannel{Tuple{S,Reversal{A}},Tuple{A}}

sensed(::FlowChannel{S,A}) where {S,A} = typeflatten(S)
actuated(::FlowChannel{S,A}) where {S,A} = typeflatten(A)

function typeflatten(s::DataType)
    map(fieldtypes(s)) do el
        el <: Tuple && return testsensed(el)
        return (el,)
    end |> Iterators.flatten |> collect |> unique!
end



abstract type Geometry end
struct NoGeometry{C} <: Geometry
    capacitance::C
end
capacitance(g::NoGeometry) = g.capacitance

abstract type PlasticityRule{S} end
struct PlasticisedChannel{S,S2,A} <: FlowChannel{Tuple{S, S2},A}
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