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
defaults(b::BasicNeuron) = merge(
    b.defaults,
    defaults(b.geometry)
)

# BasicNeuron(name::Symbol, dynamics::Dict{Q1, F}, defaults::Dict{Q2, N}, owned::Vector{C1}) where {Q1 <: Quantity, Q2 <: Quantity, F<:Function, N<:Number, C1 <: Component} = BasicNeuron{G, C1, C1, F, N, S}(name, dynamics, defaults, owned, nothing)

untagged_internal_variables(n::Neuron) = filter(x -> typeof(x) <: UntrackedQuantity, (keys(defaults(n))..., keys(dynamics(n)...) )  )

has_dynamics(n::Neuron, q::Quantity) = (q ∈ keys(dynamics(n))) ? (return true) : (return false)
has_dynamics(n::Neuron, o::OrderRelation{Q}) where Q = has_dynamics(n, Q())

liu_channels = [Liu.Na(3.0), Liu.CaS(4.0), Liu.CaT(6.0), Liu.KCa(14.0), Liu.Kdr(12.0), Liu.H(2.0), Liu.Leak(3.0)]
export liu_channels


"""
Summary:

1. build all variables required by...
- the geometry
- the subject (LHS) of the dynamical equations
- the defaults

Time dependent variable or parameter? Former if the subject of dynamical equations. Latter otherwise. Reasoning: a variable can't be dynamic if it's not the LHS of one of the dynamical equations.


2. load owned systems (both sensed and actuated by neuron being built.). Build their ODESystems.

3. hook sensed(::owned_systems) to corresponding variable of the neuron, if the owned system actually has that variable

4. sensor_defaults: for each owned system, for each sensed variable of the owned system hook the default specified by the neuron to the sensed variable.
(So the neuron requires a default value for each sensed variable of each owned system)

5. 


"""

function (b::BasicNeuron)()
   

    internal_vars = (build_vars(geometry(b)), build_vars(b |> dynamics, b), build_vars(b |> defaults, b)) |> Iterators.flatten |> Set
    # (build_vars(geometry(b))..., build_vars(b |> dynamics, b)..., build_vars(b |> defaults, b)... )
    
    all_owned = vcat(b.owned, b.preceding)
    owned_systems = [el(b) for el in all_owned]
    # preceding_systems = [el(b) for el in b.preceding]

    dynamic_equations = [f(b, owned_systems, internal_vars...) for f in values(b |> dynamics)]


    #  only hook dynamic variables. modelingtoolkit screws up for hooking parameters. hence the ∩keys(dynamics(b))
    sensor_hooks = map(all_owned, owned_systems) do channel, channel_sys 
                       (find_from(quantity, internal_vars...) ~ getproperty(channel_sys, quantity |> shorthand_name) for quantity in (sensed(channel) ∩ keys(dynamics(b))) if hasproperty(channel_sys, quantity |> shorthand_name))
        end |> Iterators.flatten |> collect


    sensor_defaults = mapreduce(merge, all_owned, owned_systems) do channel, channel_sys
        Dict(getproperty(channel_sys, quantity |> shorthand_name) => defaults(b)[quantity] for quantity in sensed(channel) if (hasproperty(channel_sys, quantity |> shorthand_name) && haskey(defaults(b), quantity)))
    end
    
    _defaults = Dict(find_from(el, internal_vars...) => defaults(b)[el] for el in (keys(defaults(b))))

    # other_defaults = ?
    return ODESystem(
        vcat(dynamic_equations, sensor_hooks),
        t;systems = owned_systems, defaults = merge(_defaults, sensor_defaults), name=b.name)
end

