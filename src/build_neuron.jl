function (ch::IonChannel)()
    @variables V(t) Ca(t)
    eqs, states, parameters, current, defaultmap = channel_dynamics(ch, V, Ca)
    sys = ODESystem(eqs, t, [V, Ca, states...], [parameters...]; observed=current,
        defaults=defaultmap, name=get_name(ch))
    return ComponentSystem(ch, sys)
end

function (ch::RegIonChannel)()
    @variables V(t) Ca(t)
    eqs, states, parameters, current, defaultmap = channel_dynamics(ch.ionch, V, Ca)
    sys = ODESystem(eqs, t, [V, states...], [parameters...]; observed=current,
        defaults=defaultmap, name=get_name(ch.ionch))
    gparam = getproperty(sys, get_g(ch.ionch); namespace=false)
    RHS = current[1].rhs
    LHS = current[1].lhs
    curr_idx = findfirst(x -> x == current[1], equations(sys))

    eqs, states, parameters, defaultmap = regul_dyn(ch, Ca)
    gstate = states[1]
    @named regul_sys = ODESystem(eqs, t, [Ca, states...], [parameters...]; defaults=defaultmap)

    equations(sys)[curr_idx] = LHS ~ gstate / gparam * RHS #take g(t) in current eq.
    sys = alias_elimination(extend(sys, regul_sys))
    return ComponentSystem(ch, sys)
end

function (syn::Synapse)(; name=get_name(syn))
    @variables Vpre(t) Vpost(t)
    eqs, states, parameters, current, defaultmap = channel_dynamics(syn, Vpre, Vpost)
    sys = ODESystem(eqs, t, [Vpre, Vpost, states...], [parameters...]; observed=current,
        name=name, defaults=defaultmap)
    return ComponentSystem(syn, sys)
end



function build_neuron(comp::Soma, channels::Vector{IonChannel}, hooks::Integer, name::Symbol)

    neuron_parameters = @parameters area Cₘ τCa Ca∞ Iapp
    neuron_states = @variables V(t) Ca(t)
    syns = [@variables $el(t) for el in [Symbol(:Isyn, i) for i = 1:hooks]]
    my_sum(syns) = hooks == 0 ? sum(Num.(syns)) : sum(reduce(vcat, syns))

    channel_systems = [ch() for ch in channels]

    summed_membrane_currents = sum(ionic_current(cs.c, cs.sys) for cs in channel_systems)
    summed_calcium_flux = sum(calcium_current(cs.c, cs.sys) for cs in channel_systems)
    connections = cat(
        [voltage_hook(V, cs) for cs in channel_systems],
        [calcium_hook(Ca, cs) for cs in channel_systems],
        dims=1
    )
    f = 14.96 # uM*nF/nA 
    #ionic currents are already carrying negative sign 
    diffeqs = cat(
        [D(V) ~ (1 / Cₘ) * (summed_membrane_currents + my_sum(syns) + Iapp)],
        [D(Ca) ~ (1 / τCa) * (-Ca + Ca∞ + (f * area * summed_calcium_flux / Cₘ))],
        dims=1
    )

    eqs = filter(x -> !(x == (Num(0) ~ Num(0))), cat(connections, diffeqs; dims=1))

    channel_defs = Dict(
        getproperty(cs.sys, el) => comp.parameters[el]
        for cs in channel_systems for el in external_params(cs.c) if (haskey(comp.parameters, el) && hasproperty(cs.sys, el))
    )

    soma_defs = Dict(
        el => comp.parameters[ModelingToolkit.tosymbol(el)] for el in neuron_parameters
    )

    state_defs = Dict(V => comp.initial_states[:V], Ca => comp.initial_states[:Ca])
    neur = ODESystem(eqs, t;
        systems=[cs.sys for cs in channel_systems],
        defaults=merge(soma_defs, channel_defs, state_defs),
        name=name
    )
    if hooks == 0
        neur = structural_simplify(neur)
    end
    return neur
end


"""
Advantage of previous method:
- never had to rebuild build_neuron, as there were the right amount of hooks
- but build_neuron(AB) is 0.014 seconds

Potential new method doesn't work
- make pre and post in add_connection each time, by build_neuron(pre_neuron)
- but need to mutably change equations from old grouping each time

1. don't add any systems when you initially build group
"""

function build_group(neurons::Vector{N}; name=:group) where {N<:Neuron}
    eqs = Array{Equation,1}(undef, 0)
    ODESystem(eqs, t, [], []; name=name, systems=[el.ODESystem for el in neurons])
end


"""
Issue: need to do all connections in one ḡChol

"""

function add_connection(group, pre, post, syn::EmptyConnection; kwargs...)
    return group
end

function add_connection(group, pre_n, post_n, syn::Synapse; name=ModelingToolkit.getname(group), i=1)
    pre = pre_n.ODESystem
    post = post_n.ODESystem
    prename, postname = ModelingToolkit.getname.([pre, post])
    synapse_sys = syn(; name=Symbol(prename, :to, postname, get_name(syn)))

    oldeqs = ModelingToolkit.get_eqs(group)
    neweqs = [
        get_pre(syn, pre) ~ my_pre(synapse_sys),
        get_post(syn, post) ~ my_post(synapse_sys),
        syn_current(syn, synapse_sys.sys) ~ getproperty(post, Symbol(post_connector(syn), i))
    ]

    eqs = cat(oldeqs, neweqs, dims=1)
    systems = cat(ModelingToolkit.get_systems(group), synapse_sys.sys, dims=1)

    return connected = ODESystem(eqs, t, [], [];
        name=name,
        systems=systems)
end

