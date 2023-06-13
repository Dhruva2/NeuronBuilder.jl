
"""
NB in further instantiations, the following might be our friends:
ModelingToolkit.connect 
ModelingTOolkti.unpack



requirements:


Vpre and Vpost. how? 

Need a naming and type convention:

Pre{Quantity}
-> sys.xxxpre

Post{Qyantity}
-> sys.xxxpost


for el in Pre
    sys.el ~ pre_sys.el # sensed parts
end

for el in Post
    sys.el ~ post_sys.el # sensed parts #NOT INCLUDING ACTUATIONS
end

actuations are hooked in post_neuron building script. how to distoinguish Current{Voltage} and Current{Post{Voltage}}? currently eg Cdv/dt dynamics only looks for Current{Voltage}


Option 1: get rid of Post. Only have Pre. since the synapse anyway beongs to the actuated channel.
Much easier than some weird stripping of Post at particular points.
Any issues? Can't think of any



1. 
"""

"""
Basic synapse {Cholinergic}

sensed(Previous, s::Synapse) = (Voltage(), )

ISSUE: how do I determine if sensed quantities are dynaimc? 
CURRENT WORKAROUND: set this manually for now. BasicSynapses only sense voltage, and this is dynamic.

"""

function build_vars_from_owner(s::DirectedChannel, post::Component, which::Function)
    return Iterators.map(which(s)) do quantity
        _name = shorthand_name(quantity)
        if typeof(quantity) <: Previous || has_dynamics(post, quantity)
            el, = @variables $_name(t)
        else
            el, = @parameters $_name
        end
        return el
    end
end



struct BasicSynapse{Q1<:Quantity,Q2<:Quantity,Q3<:Quantity,F<:Function,N} <: DirectedChannel
    name::Symbol
    type::Q1
    dynamics::Dict{Q2, F}
    defaults::Dict{Q3, N}
end


pre_sensed(b::BasicSynapse) = Set((Previous{Voltage}(),))
post_sensed(b::BasicSynapse) = Set((Voltage(),))
sensed(b::BasicSynapse) = pre_sensed(b) ∪ post_sensed(b)
actuated(b::BasicSynapse) = Set((Current{Voltage}(),))
tagged_internal_variables(b::BasicSynapse) = Set((Conductance{b.type |> typeof}(), Reversal{b.type |> typeof}()))





function (ch::BasicSynapse)(owner::Component)
    vars = get_all_vars(ch, owner)
    eqs = [v(ch, vars...) for v in values(ch |> dynamics)] #|> Iterators.flatten |> collect
    _defaults = Dict(find_from(el, vars...) => defaults(ch)[el] for el in (keys(defaults(ch))) if typeof(defaults(ch)[el]) <: Number)
    calculated_defaults = Dict(find_from(el, vars...) => defaults(ch)[el](vars...) for el in (keys(defaults(ch))) if typeof(defaults(ch)[el]) <: Function)

    return ODESystem(eqs, t, defaults=merge(_defaults, calculated_defaults), name=ch.name)
end






# _dynamics = Dict(
#     UntrackedQuantity(:s) => Liu.synaptic_channel_dynamics[:s],
#     Current{Voltage}() => BasicComponents.basic_mh_current(Choline(), :s => 1, :h => 0; output = Voltage())
# )

# _defaults = Dict(
#     UntrackedQuantity(:k₋) => 0.025,
#     UntrackedQuantity(:Vth) => -35.0,
#     UntrackedQuantity(:δ) => 5.0,
#     UntrackedQuantity(:s) => 0.0,
#     Reversal{Choline}() => -70.0
# )

# chol = BasicSynapse(
#     :chol,
#     Choline(),
#     _dynamics,
#     _defaults
# )





# eqs = [D(s) ~ (1 / τs(ch, Vpre)) * (s̄(ch, Vpre) - s),
#     IChol ~ -ḡChol * s * (Vpost - ch.Eₛ)]
# current = [eqs[2]]
