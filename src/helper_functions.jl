#return unparameterised type as symbol
function get_name(ch::Component)
    Base.typename(ch |> typeof).name |> Symbol
end

function get_name(ch)
    Base.typename(ch).name |> Symbol
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
shorthand_name(::Type{PseudoIon}) = :Leak
shorthand_name(::Type{Reversal{T}}) where {T} = Symbol(:E, shorthand_name(T))
shorthand_name(::Type{Current{T}}) where {T} = Symbol(:I, shorthand_name(T))
shorthand_name(::Type{Conductance{T}}) where {T} = Symbol(:g, shorthand_name(T))
shorthand_name(::Type{mRNA{T}}) where {T} = Symbol(:mRNA_, shorthand_name(T))

shorthand_name(x::Type{Tuple{T,R}}) where {T,R} = shorthand_name.(x.types)


sensed(::FlowChannel{S,A}) where {S,A} = typeflatten(S)
actuated(::FlowChannel{S,A}) where {S,A} = typeflatten(A)





function typeflatten(s::DataType)
    map(fieldtypes(s)) do el
        el <: Tuple && return typeflatten(el)
        return (el,)
    end |> Iterators.flatten |> collect |> unique!
end

voltage(el) = (Voltage,)
ionic(x) = x <: Ion

currents(i::Channel) =
    map(filter(ionic, actuated(i))) do thing
        Current{thing}
    end

sensedvars(i::Channel) =
    map(sensed(i)) do thing
        Symbol(thing |> shorthand_name)
    end

mrna(i::Channel) =
    map(filter(ionic, actuated(i))) do thing
        mRNA{thing}
    end

reversals(i::Channel) =
    map(filter(ionic, actuated(i))) do thing
        Reversal{thing}
    end

conductances(i::Channel) =
    map(filter(ionic, actuated(i))) do thing
        Conductance{thing}
    end

sensed_ions(i::Channel) = filter(ionic, sensed(i))





pre_sensed(::Synapse{Spre,Spost,A}) where {Spre,Spost,A} = typeflatten(Spre)
pre_sensed(::Synapse{Spre,Spost,A}) where {Spre,Spost,A} = typeflatten(Spost)
actuated(::Synapse{Spre,Spost,A}) where {Spre,Spost,A} = typeflatten(A)












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


instantiate_variables(v::Vector{Symbol}) =
    map(v) do f
        @variables $f(t)
    end |> Iterators.flatten |> collect


instantiate_parameters(v::Vector{Symbol}) =
    map(v) do f
        @variables $f
    end |> Iterators.flatten |> collect


get_actuator(c::FlowChannel, sys::ODESystem, v::Type{Voltage}) =
    sum(currents(c)) do I
        ModelingToolkit.getvar(sys, shorthand_name(I); namespace=true)
    end

function get_actuator(c::FlowChannel, sys::ODESystem, f::DataType)
    indx = findfirst(x -> x == f, actuated(c))
    isnothing(indx) && return Num(0.0)
    return ModelingToolkit.getvar(sys, shorthand_name.(currents(c))[indx]; namespace=true)
end

function get_sensor(c::FlowChannel, sys::ODESystem, f::DataType)
    if f ∈ sensed(c)
        return ModelingToolkit.getvar(sys, shorthand_name(f); namespace=true)
    end
end
get_sensor(c::FlowChannel, sys::ODESystem, f::Type{Voltage}) =
    ModelingToolkit.getvar(sys, shorthand_name(f); namespace=true)

get_parameters(::SpeciesDynamics) = nothing
get_states(::SpeciesDynamics) = nothing
default_params(::SpeciesDynamics, a, b, c) = Dict{Num,Float64}()
default_states(::SpeciesDynamics, a, b, c) = Dict{Num,Float64}()


function get_from(d::Dict{DataType,SpeciesDynamics}, func)
    muddled = d |> values .|> func |> x -> filter(y -> !isnothing(y), x)
    isempty(muddled) && return Vector{Num}[]
    return reduce(vcat, muddled)::Vector{Num}
end

function vardivide(v::Num...)
    states = [filter(!ModelingToolkit.isparameter, v)...]
    params = [filter(ModelingToolkit.isparameter, v)...]
    return states, params
end
