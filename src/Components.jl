abstract type AbstractComponent end 
abstract type Component <: AbstractComponent end

function get_name(ch::Component)
    Base.typename(ch |> typeof).name |> Symbol
end


abstract type Compartment <: Component end
abstract type BasicChannel <: Component end
abstract type DirectedChannel <: BasicChannel end
struct EmptySynapse <: DirectedChannel end 



function is_dynamic(b::BasicChannel, q::Quantity)
    (q ∈ keys(b.dynamics)) ? (return true) : (return false)
end

is_dynamic(q::Quantity, n::Neuron) = (q ∈ keys(dynamics(n)))

abstract type Neuron <: Compartment end
abstract type PlasticityRule end
abstract type Geometry <: AbstractComponent end


struct NoGeometry{Q<:Quantity, N<:Number} <: Geometry
    defaults::Dict{Q, N}
end

dynamics(::NoGeometry) = Dict()
defaults(n::NoGeometry) = n.defaults

parameters(g::Geometry) = Dict(el => defaults(g)[el] for el in setdiff(keys(defaults(g)), keys(dynamics(g))))
is_dynamic(q::Quantity, g::Geometry) = (q ∈ keys(dynamics(g)))
build_vars(g::Geometry) = (build_vars(merge(dynamics(g), parameters(g)), g)...,)





### if you want to build a channel without these fields then overwrite the corresponding functions 
dynamics(b::BasicChannel) = b.dynamics
name(b::BasicChannel) = b.name
defaults(b::BasicChannel) = b.defaults 


function build_vars_from_owner(comp::Component, owner::Component, which::Function)
    return map(which(comp)) do quantity
        _name = shorthand_name(quantity)
        if has_dynamics(owner, quantity)
            @variables $_name(t)
        else
            @parameters $_name
        end
    end |> Iterators.flatten |> collect
end

build_vars(d::Dict, c::AbstractComponent) = map(keys(d) |> collect) do el
        _name = shorthand_name(el)
        if is_dynamic(el, c)
            @variables $_name(t)
        else
            @parameters $_name
        end
end |> Iterators.flatten |> collect


function build_vars(comp::Component, which::Function)
    return map(which(comp)) do quantity
               _name = shorthand_name(quantity)
               if is_dynamic(comp, quantity)
                   @variables $_name(t)
               else
                   @parameters $_name
               end
           end |> Iterators.flatten |> collect
end

function build_vars(q::Quantity...)
     _name = shorthand_name(quantity)
end

function build_defaults(c::Component, sys::ModelingToolkit.AbstractSystem)
    return Dict(
                getproperty(sys, shorthand_name(key)) => val 
                for (key, val) in defaults(c)
                    )
end


function find_from(q::Quantity, vars...)
    _name = shorthand_name(q)
    for el in vars
        if _name == ModelingToolkit.tosymbol(el, escape=false)
            return el
        end
    end
end

function (ch::BasicChannel)(owner::Component)
    vars = get_all_vars(ch, owner)
    eqs = [v(ch, vars...) for v in values(ch |> dynamics)] |> Iterators.flatten |> collect

    _defaults = Dict(find_from(el, vars...) => defaults(ch)[el] for el in (keys(defaults(ch))) if typeof(defaults(ch)[el])<:Number)

    calculated_defaults = Dict(find_from(el, vars...) => defaults(ch)[el](vars...) for el in (keys(defaults(ch))) if typeof(defaults(ch)[el])<:Function)

    return  ODESystem(eqs, t, defaults=merge(_defaults, calculated_defaults), name=ch.name)
end


function get_all_vars(ch::Component, owner::Component)
    vars = (
        build_vars_from_owner(ch, owner, sensed)...,
        build_vars(ch, tagged_internal_variables)...,
        build_vars(ch, untagged_internal_variables)...,
        build_vars(ch, actuated)...
    )
end


function find_from(q::Vector{Q}, vars...) where Q <: Quantity
    return [find_from(el, vars...) for el in q]
end


function calc_defaults_from(f::Function, q::Vector{Q}) where Q<:Quantity
    
    return function _f(vars...)
        inputs = find_from(q, vars...)
        return f(inputs...)
    end
end
export calc_defaults_from


