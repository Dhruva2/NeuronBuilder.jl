function (ch::IonChannel)()
    @variables V(t) Ca(t)
    eqs, states, parameters, current, defaultmap = channel_dynamics(ch, V, Ca)
    sys = ODESystem(eqs, t, [V, Ca, states...], [parameters...];
        observed = current, defaults = defaultmap, name = get_name(ch))
    return ComponentSystem(ch, sys)
end

function (ch::RegIonChannel)()
    @variables V(t) Ca(t)
    eqs, states, parameters, current, defaultmap = channel_dynamics(ch.ionch, V, Ca)
    sys = ODESystem(eqs, t, [V, states...], [parameters...];
        observed = current, defaults = defaultmap, name = get_name(ch.ionch))
    RHS = current[1].rhs
    LHS = current[1].lhs
    curr_idx = findfirst(x -> x == current[1], equations(sys))

    eqs, states, parameters, defaultmap = regul_dyn(ch, Ca)
    @named regul_sys = ODESystem(eqs, t, [Ca, states...], [parameters...]; defaults = defaultmap)

    gstate = regul_sys.states[2]
    gparam = sys.ps[1]
    equations(sys)[curr_idx] = LHS ~ gstate / gparam * RHS #take g(t) in current eq.
    sys = alias_elimination(extend(sys, regul_sys))
    return ComponentSystem(ch.ionch, sys)
end

function (syn::Synapse)(; name = get_name(syn))
    @variables Vpre(t) Vpost(t)
    eqs, states, parameters, current, defaultmap = channel_dynamics(syn, Vpre, Vpost)
    sys = ODESystem(eqs, t, [Vpre, Vpost, states...], [parameters...]; observed = current,
        name = name, defaults = defaultmap)
    return ComponentSystem(syn, sys)
end


function build_neuron(comp, channels; name = :unidentified_neuron, reg = false)

    neuron_parameters = @parameters Cₘ τCa Ca∞ Iapp
    neuron_states = @variables V(t) Ca(t)
    syns = [@variables $el(t) for el in [Symbol(:Isyn, i) for i = 1:comp.hooks]]
    my_sum(syns) = comp.hooks == 0 ? sum(Num.(syns)) : sum(reduce(vcat, syns))

    if reg
        channels = [isLeak(ch) ? ch : RegIon(ch, getfield(ch, get_g(ch)), Symbol(:τ, get_name(ch))) for ch in channels]
    end
    channel_systems = [ch() for ch in channels]

    summed_membrane_currents = sum(ionic_current(cs.c, cs.sys) for cs in channel_systems)

    summed_calcium_flux = sum(calcium_current(cs.c, cs.sys) for cs in channel_systems)

    connections = cat(
        [voltage_hook(V, cs) for cs in channel_systems],
        [calcium_hook(Ca, cs) for cs in channel_systems],
        dims = 1
    )
    f = 0.94 #ionic currents are already carrying negative sign
    diffeqs = cat(
        [D(V) ~ (summed_membrane_currents + my_sum(syns) + Iapp) / Cₘ],
        [D(Ca) ~ (1 / τCa) * (-Ca + Ca∞ + (f * summed_calcium_flux))],
        dims = 1
    )

    eqs = cat(connections, diffeqs; dims = 1)
    eqs = filter(x -> !(x == (Num(0) ~ Num(0))), eqs)

    channel_defs = Dict{Num,Float64}()
    for cs in channel_systems
        param_names = filter(x -> x !== nothing, external_params(cs.c))
        for el in param_names
            if haskey(comp.parameters, el) && hasproperty(cs.sys, el)
                channel_defs[getproperty(cs.sys, el)] = comp.parameters[el]
            end
        end
    end

    soma_defs = Dict{Num,Float64}()
    for el in neuron_parameters
        soma_defs[el] = comp.parameters[ModelingToolkit.tosymbol(el)]
    end
    state_defs = Dict(V => comp.initial_states[:V], Ca => comp.initial_states[:Ca])
    defaults = merge(soma_defs, channel_defs, state_defs)
    neur = ODESystem(eqs, t;
        systems = [cs.sys for cs in channel_systems],
        defaults = defaults,
        name = name
    )

    return neur
end


function build_group(arr_of_neurons; name = :group)
    eqs = Array{Equation,1}(undef, 0)
    ODESystem(eqs, t, [], []; name = name, systems = arr_of_neurons)
end


"""
Issue: need to do all connections in one ḡChol

"""
function add_connection(group, pre, post, syn::Synapse; name = ModelingToolkit.getname(group), i = 1)
    prename, postname = ModelingToolkit.getname.([pre, post])
    synapse_sys = syn(; name = Symbol(prename, :to, postname, get_name(syn)))

    oldeqs = ModelingToolkit.get_eqs(group)
    neweqs = [
        get_pre(syn, pre) ~ my_pre(synapse_sys),
        get_post(syn, post) ~ my_post(synapse_sys),
        syn_current(syn, synapse_sys.sys) ~ getproperty(post, Symbol(post_connector(syn), i))
    ]

    eqs = cat(oldeqs, neweqs, dims = 1)
    systems = cat(ModelingToolkit.get_systems(group), synapse_sys.sys, dims = 1)

    return connected = ODESystem(eqs, t, [], [];
        name = name,
        systems = systems)
end

