abstract type Component end

function get_name(ch::Component)
    Base.typename(ch |> typeof).name |> Symbol
end




abstract type Compartment <: Component end
abstract type FlowChannel{Sensors<:Tuple,Actuators<:Tuple} <: Component end
abstract type Neuron <: Compartment end


abstract type Synapse <: Component end
struct EmptyConnection <: Synapse end

abstract type Geometry end


abstract type PlasticityRule{S} end

struct PlasticisedChannel{S,S2,A,C<:FlowChannel{S,A},P<:PlasticityRule{S2}} <: FlowChannel{Tuple{S,S2},A}
    channel::C
    modification::P
end

function get_name(p::PlasticisedChannel)
    Symbol(
        get_name(p.channel),
        :_with_,
        get_name(p.mutation)
    )
end

