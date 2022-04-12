abstract type Species end
abstract type Ion <: Species end
abstract type AbstractIon <: Ion end

"""
Not including voltage in tagging ion channels. Any flow channel. Flows are all charged particles. So they will always include voltage. Not going to consider flows of e.g. proteins
"""
struct Voltage <: Species end
struct Sodium <: Ion end
struct Potassium <: Ion end
struct Calcium <: Ion end
struct Proton <: Ion end
struct Leak <: AbstractIon end
abstract type Reversal{<:Ion} end


function merge_types(t1::Type, t2::Type)
    ps = (type_ps(t1)..., type_ps(t2)...) |> unique
    return Tuple{ps...}
end

function type_ps(t::Type)
    isempty(t.parameters) && return (t,)
    return t.parameters |> Tuple
end

abstract type Top{A<:Tuple,B<:Tuple} end

struct Example{T} <: Top{Tuple{Int64},Tuple{DataType}}
    a::T
end

struct __MutatedExample{A,A2<:Tuple,B} <: Top{merge_types(Type{A}, Type{A2}),B}
    e::Top{A,B} # in particular, e::Example{T}
    mutation::A2
end

function MutatedExample(e::Top{A,B}, m::A2) where {A,B,A2}
    MutatedExample{A,T,A2,B}(e, m)
end

struct FlowChannel
    sensed
    actuated
end


cu