"""

abstract type FlowChannel{Sensors<:Tuple,Actuators<:Tuple} <: Component end

"""




# channels always sense the reversal of the currents they actuate
FlowChannel(T) = FlowChannel{Tuple{Reversal{T}},Tuple{T}}
FlowChannel(S, A) = FlowChannel{Tuple{S,Reversal{A}},Tuple{A}}
FlowChannel() = FlowChannel{Tuple{Nothing},Tuple{Nothing}}


sensed(::FlowChannel{S,A}) where {S,A} = typeflatten(S)
actuated(::FlowChannel{S,A}) where {S,A} = typeflatten(A)
sensed(::FlowChannel{Tuple{Nothing},A}) where {A} = Vector{DataType}()
actuated(::FlowChannel{S,Tuple{Nothing}}) where {S} = Vector{DataType}()

function typeflatten(s::DataType)
    map(fieldtypes(s)) do el
        el <: Tuple && return typeflatten(el)
        return (el,)
    end |> Iterators.flatten |> collect |> unique!
end

voltage(el) = (Voltage,)
ionic(x) = x <: Ion

currents(i::FlowChannel) =
    map(filter(ionic, actuated(i))) do thing
        Current{thing}
    end

sensedvars(i::FlowChannel) =
    map(sensed(i)) do thing
        Symbol(thing |> shorthand_name)
    end

mrna(i::FlowChannel) =
    map(filter(ionic, actuated(i))) do thing
        mRNA{thing}
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


instantiate_variables(v::Vector{Symbol}) =
    map(v) do f
        @variables $f(t)
    end |> Iterators.flatten |> collect

instantiate_variables(s::Symbol) = @variables $s(t)

instantiate_parameters(v::Vector{Symbol}) =
    map(v) do f
        @variables $f
    end |> Iterators.flatten |> collect

instantiate_parameters(s::Symbol) = @parameters $f


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