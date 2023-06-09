"""
tagged internal variables:
    everything sensed
    NOT everything actuated
    every key in n.dynamics

no sensed / actuated

untagged internal variables:
capacitance, anything in geometry 


has_dynamics:
every key in n.dynamics
anything else? don't think so
"""


struct BasicNeuron{G<:Geometry,Q1 <: Quantity, F<:Function, Q2 <: Quantity, N<:Number, C1<:BasicChannel,C2<:DirectedChannel} <: Neuron
    name::Symbol
    geometry::G
    dynamics::Dict{Q1,F}
    defaults::Dict{Q2, N}
    owned:: Vector{C1}
    preceding::Vector{C2}
end

geometry(b::BasicNeuron) = b.geometry
dynamics(b::BasicNeuron) = b.dynamics
defaults(b::BasicNeuron) = b.defaults

# BasicNeuron(name::Symbol, dynamics::Dict{Q1, F}, defaults::Dict{Q2, N}, owned::Vector{C1}) where {Q1 <: Quantity, Q2 <: Quantity, F<:Function, N<:Number, C1 <: Component} = BasicNeuron{G, C1, C1, F, N, S}(name, dynamics, defaults, owned, nothing)

untagged_internal_variables(n::Neuron) = filter(x -> typeof(x) <: UntrackedQuantity, (keys(defaults(n))..., keys(dynamics(n)...) )  )

has_dynamics(n::Neuron, q::Quantity) = (q ∈ keys(dynamics(n))) ? (return true) : (return false)

liu_channels = [Liu.Na(3.0), Liu.CaS(4.0), Liu.CaT(6.0), Liu.KCa(14.0), Liu.Kdr(12.0), Liu.H(2.0), Liu.Leak(3.0)]
export liu_channels




export geometry, dynamics, defaults, build_defaults

"""
what variables are missing from map?

Those that are sensed by the owned, and defaulted from the owner instead of the owned

we SHOULD NOT instantiate the internal variables ENa etc. ie the parameters used by the owned

some of the above might be needed by dynamic equations. 
so we should NOT include them in the overall sys equations through sensor hooks. 


"""

function (b::BasicNeuron)()
   

    internal_vars = (build_vars(geometry(b))..., build_vars(b |> dynamics, b)..., build_vars(b |> defaults, b)... )
    owned_systems = [el(b) for el in b.owned]
    dynamic_equations = [f(b, owned_systems, internal_vars...) for f in values(b |> dynamics)]


    #  only hook dynamic variables. modelingtoolkit screws up for hooking parameters. hence the ∩keys(dynamics(b))
    sensor_hooks = map(b.owned, owned_systems) do channel, channel_sys 
                       [find_from(quantity, internal_vars...) ~ getproperty(channel_sys, quantity |> shorthand_name) for quantity in (sensed(channel) ∩ keys(dynamics(b))) if hasproperty(channel_sys, quantity |> shorthand_name)]
        end |> Iterators.flatten |> collect


    sensor_defaults = mapreduce(merge, b.owned, owned_systems) do channel, channel_sys
        Dict(getproperty(channel_sys, quantity |> shorthand_name) => defaults(b)[quantity] for quantity in sensed(channel) if hasproperty(channel_sys, quantity |> shorthand_name))
    end
    
    _defaults = Dict(find_from(el, internal_vars...) => defaults(b)[el] for el in (keys(defaults(b))))
    _geom_defaults = Dict(find_from(el, internal_vars...) => defaults(b |>geometry)[el] for el in (keys(defaults(b |>geometry))))

    # other_defaults = ?

    return ODESystem(
        vcat(dynamic_equations, sensor_hooks),
        t;systems = owned_systems, defaults = merge(_defaults, _geom_defaults, sensor_defaults), name=b.name)
end

