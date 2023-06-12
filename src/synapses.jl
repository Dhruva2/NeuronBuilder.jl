"""
Previous system:

- Currents totted up manually in neuron()

Hooking:
-add_connection() did pre ~ pre and post ~ post only for voltage


"""


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
    sys.el ~ pre_sys.el # sensed parts #NOT INCLUDING ACTUATIONS
end

actuations are hooked in post_neuron building script. 

1. 
"""

"""
Basic synapse {Cholinergic}

sensed(Previous, s::Synapse) = (Voltage(), )

"""



struct BasicSynapse{Q1<:Quantity,Q2<:Quantity,Q3<:Quantity,F<:Function,N} <: DirectedChannel
    name::Symbol
    type::Q1
    dynamics::Dict{Q2, F}
    defaults::Dict{Q3, N}
end


pre_sensed(b::BasicSynapse) = Set(Voltage(),)
post_sensed(b::BasicSynapse) = Set(Voltage(),)
sensed(b::BasicSynapse) = pre_sensed(b) ∪ post_sensed(b)
# sensed(b::BasicSynapse) = ((pre_sensed(b) ∪ post_sensed(b))...,)
actuated(b::BasicSynapse) = Set((Current{Voltage}(),))
tagged_internal_variables(b::BasicSynapse) = Set(Conductance{b.type})







function (ch::BasicSynapse)(owner::Component)
    vars = get_all_vars(ch, owner)
    eqs = [v(ch, vars...) for v in values(ch |> dynamics)] |> Iterators.flatten |> collect
    _defaults = Dict(find_from(el, vars...) => defaults(ch)[el] for el in (keys(defaults(ch))) if typeof(defaults(ch)[el]) <: Number)
    calculated_defaults = Dict(find_from(el, vars...) => defaults(ch)[el](vars...) for el in (keys(defaults(ch))) if typeof(defaults(ch)[el]) <: Function)

    return ODESystem(eqs, t, defaults=merge(_defaults, calculated_defaults), name=ch.name)
end



"""
dynamics: s, Ichol
parameters: k₋, Vth, δ
"""

_defaults = Dict(
    UntrackedQuantity(:k₋)  => 0.025,
    UntrackedQuantity(:Vth) => -35.0,
    UntrackedQuantity(:δ)   => 5.0,
    UntrackedQuantity(:s̄)   => 0.0,
    UntrackedQuantity(:Eₛ)  => -70.0
)

_dynamics = Dict(

)

UntrackedQuantity(:s̄) => 1.0 / (1.0 + exp((syn.Vth - Vpre) / syn.δ))

τs(syn::Glut, Vpre) = (1.0 - s̄(syn, Vpre)) / syn.k₋

function Cholinergic_dynamics(c::Component, vars...)
    τ(s) = 
        s̄, k₋, Vth, δ = find_from(c, [UntrackedQuantity(:s̄), UntrackedQuantity(:k₋)], vars...)
    return [D(s) ~ ...]
end




dynamics = Dict(

)

chol = BasicSynapse(
    :chol,
    Choline(),
    s,
    b
)

