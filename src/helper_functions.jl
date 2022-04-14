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

sensed(::FlowChannel{S,A}) where {S<:Tuple,A} = fieldtypes(S)
actuated(::FlowChannel{S,A}) where {S,A<:Tuple} = fieldtypes(A)

sensed(::FlowChannel{S,A}) where {S,A} = (S,)
actuated(::FlowChannel{S,A}) where {S,A} = (A,)


sensedvars(i::FlowChannel) =
    map(sensed(i)) do thing
        Symbol(thing |> shorthand_name)
    end

reversals(i::FlowChannel) =
    map(actuated(i)) do thing
        Symbol(:E, shorthand_name(thing))
    end

currents(i::FlowChannel) =
    map(actuated(i)) do thing
        Symbol(:I, shorthand_name(thing))
    end

conductance(i::FlowChannel) =
    map(actuated(i)) do thing
        Symbol(:g, shorthand_name(thing))
    end

instantiate_variables(c::Component, args...) =
    map(args) do el
        map(el(c)) do f
            @variables $f(t)
        end |> Iterators.flatten |> collect
    end |> Iterators.flatten |> collect

instantiate_parameters(c::Component, args...) =
    map(args) do el
        map(el(c)) do f
            @parameters $f
        end |> Iterators.flatten |> collect
    end |> Iterators.flatten |> collect


instantiate_variables(v::Vector{Symbol}) =
    map(v) do f
        @variables $f(t)
    end |> Iterators.flatten |> collect

get_actuator(c::ComponentSystem{C,S}, v::Type{Voltage}) where {C<:FlowChannel,S} =
    sum(currents(c.c)) do I
        ModelingToolkit.getvar(c.sys, I; namespace=true)
    end

function get_actuator(c::ComponentSystem{C,S}, f::DataType) where {C<:FlowChannel,S}
    indx = findfirst(x -> x == f, actuated(c.c))
    isnothing(indx) && return Num(0.0)
    return ModelingToolkit.getvar(c.sys, currents(c.c)[indx]; namespace=true)
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

#isLeak(c::IonChannel) = typeof(c) in [NeuronBuilder.Liu.Leak{Float64}, NeuronBuilder.Prinz.Leak{Float64}]