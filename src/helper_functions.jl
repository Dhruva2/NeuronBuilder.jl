######### helper functions ###########

#return unparameterised type as symbol
function get_name(ch::Component)
    Base.typename(ch |> typeof).name |> Symbol
end

function get_name(p::PlasticisedChannel)
    Symbol(
        get_name(p.channel),
        :_with_,
        get_name(p.mutation)
    )
end

shorthand_name(::Type{Voltage}) = :V
shorthand_name(::Type{Sodium}) = :Na
shorthand_name(::Type{Potassium}) = :K
shorthand_name(::Type{Calcium}) = :Ca
shorthand_name(::Type{Proton}) = :H
shorthand_name(::Type{Leak}) = :Leak
shorthand_name(::Type{Reversal{T}}) where {T} = Symbol(:E, shorthand_name(T))
shorthand_name(x::Type{Tuple{T,R}}) where {T,R} = shorthand_name.(x.types)

sensed(::FlowChannel{S,A}) where {S,A} = typeflatten(S)
actuated(::FlowChannel{S,A}) where {S,A} = typeflatten(A)

function typeflatten(s::DataType)
    map(fieldtypes(s)) do el
        el <: Tuple && return typeflatten(el)
        return (el,)
    end |> Iterators.flatten |> collect |> unique!
end

sensedvars(i::FlowChannel) =
    map(sensed(i)) do thing
        Symbol(thing |> shorthand_name)
    end

reversals(i::FlowChannel) =
    map(filter(ionic, actuated(i))) do thing
        Reversal{thing}
    end

conductances(i::FlowChannel) =
    map(filter(ionic, actuated(i))) do thing
        Conductance{thing}
    end

sensed_ions(i::FlowChannel) = filter(ionic, sensed(i))

instantiate_hooks(pre::Compartment, to::Component, args::Function...) =
    map(args) do fun
        map(fun(to)) do thing
            _name = shorthand_name(thing)
            if has_dynamics(pre, thing)
                return @variables $_name(t)
            else
                return @parameters $_name
            end
        end |> Iterators.flatten |> collect
    end |> Iterators.flatten |> collect

instantiate_variables(c::Component, args::Function...) =
    map(args) do fun
        map(fun(c)) do thing
            _name = shorthand_name(thing)
            return @variables $_name(t)
        end |> Iterators.flatten |> collect
    end |> Iterators.flatten |> collect

instantiate_parameters(c::Component, args::Function...) =
    map(args) do fun
        map(fun(c)) do thing
            _name = shorthand_name(thing)
            return @parameters $_name
        end |> Iterators.flatten |> collect
    end |> Iterators.flatten |> collect


export currents, reversals, conductances, instantiate_hooks, species, sensed_ions, voltage

instantiate_variables(v::Vector{Symbol}) =
    map(v) do f
        @variables $f(t)
    end |> Iterators.flatten |> collect


instantiate_parameters(v::Vector{Symbol}) =
    map(v) do f
        @variables $f
    end |> Iterators.flatten |> collect


get_actuator(c::ComponentSystem{C,S}, v::Type{Voltage}) where {C<:FlowChannel,S} =
    sum(currents(c.c)) do I
        ModelingToolkit.getvar(c.sys, shorthand_name(I); namespace=true)
    end

function get_actuator(c::ComponentSystem{C,S}, f::DataType) where {C<:FlowChannel,S}
    indx = findfirst(x -> x == f, actuated(c.c))
    isnothing(indx) && return Num(0.0)
    return ModelingToolkit.getvar(c.sys, shorthand_name.(currents(c.c))[indx]; namespace=true)
end

function get_sensor(c::ComponentSystem{C,S}, f::DataType) where {C<:FlowChannel,S}
    if f ∈ sensed(c.c)
        return ModelingToolkit.getvar(c.sys, shorthand_name(f); namespace=true)
    end
end
get_sensor(c::ComponentSystem{C,S}, f::Type{Voltage}) where {C<:FlowChannel,S} =
    ModelingToolkit.getvar(c.sys, shorthand_name(f); namespace=true)

get_parameters(::SpeciesDynamics) = nothing
get_states(::SpeciesDynamics) = nothing
default_params(::SpeciesDynamics, a, b, c) = Dict{Num,Float64}()
default_states(::SpeciesDynamics, a, b, c) = Dict{Num,Float64}()


function get_from(d::Dict{DataType,SpeciesDynamics}, func)
    muddled = d |> values .|> func |> x -> filter(y -> !isnothing(y), x)
    isempty(muddled) && return Vector{Num}[]
    return reduce(vcat, muddled)::Vector{Num}
end

isLeak(c::FlowChannel) = typeof(c) in [NeuronBuilder.Liu.leak{Float64}, NeuronBuilder.Prinz.leak{Float64}]

function vardivide(v::Num...)
    states = [filter(!ModelingToolkit.isparameter, v)...]
    params = [filter(ModelingToolkit.isparameter, v)...]
    return states, params
end

export vardivide

