function (ch::IonChannel)()
    @variables V(t) Ca(t)
    eqs, states, parameters, current, defaultmap = channel_dynamics(ch, V, Ca)
    sys = ODESystem(eqs, t, [V, Ca, states...], [parameters...];
        observed = current, defaults = defaultmap, name = get_name(ch))
    return ComponentSystem(ch,sys)
end

function (ch::IonChannel)(reg::Bool)
    cs = ch()
    sys = cs.sys
    if reg && typeof(ch) != Leak{Float64}
        regul_sys = regul_dyn(ch, sys.Ca)
        sys = extend(sys, regul_sys)
        push!(equations(sys), sys.parameters[1] ~ regul_sys.states[1])
    end
    return ComponentSystem(ch, sys)
end

function build_channel(syn::Synapse; name = get_name(syn))
    @variables Vpre(t) Vpost(t)
    eqs, states, parameters, current, defaultmap = channel_dynamics(syn, Vpre, Vpost)
    return ODESystem(eqs, t, [Vpre, Vpost, states...], [parameters...]; observed = current,
        name = name, defaults = defaultmap)
end


# 10 nF/mm² for default specific membrane capacitance
function build_neuron(comp, channels, reg; name = :unidentified_neuron)

    neuron_parameters = @parameters Cₘ τCa Ca∞ Ca_tgt τg
    neuron_states = @variables V(t) Ca(t)
    syns = [@variables $el(t) for el in [Symbol(:Isyn, i) for i = 1:comp.hooks]]
    my_sum(syns) = comp.hooks == 0 ? sum(Num.(syns)) : sum(syns)

    channel_systems = [ch() for ch in channels]

    summed_membrane_currents = sum(ionic_current(ch, cs.sys) for (ch, cs) in zip(channels, channel_systems))

    summed_calcium_flux = sum(calcium_current(ch, cs.sys) for (ch, cs) in zip(channels, channel_systems))

    connections = cat(
        [voltage_hook(V, cs) for cs in channel_systems],
        [calcium_hook(Ca, cs) for cs in channel_systems],
        dims = 1
    )
    f = 0.094
    diffeqs = cat(
        [D(V) ~ (summed_membrane_currents + my_sum(syns)) / Cₘ],
        [D(Ca) ~ (1 / τCa) * (-Ca + Ca∞ + f * summed_calcium_flux)],
        dims = 1
    )

    eqs = cat(connections, diffeqs; dims = 1)
    eqs = filter(x -> !(x == (Num(0) ~ Num(0))), eqs)

    channel_defs = Dict{SymbolicUtils.Symbolic{Real},Float64}()
    for cs in channel_systems
        param_names = external_params(cs.c)
        if !isnothing(param_names)
            for el in param_names
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
    synapse = build_channel(syn; name = Symbol(prename, :to, postname, get_name(syn)))

    oldeqs = ModelingToolkit.get_eqs(group)
    neweqs = [
        get_pre(syn, pre) ~ my_pre(syn, synapse),
        get_post(syn, post) ~ my_post(syn, synapse),
        syn_current(syn, synapse) ~ getproperty(post, Symbol(post_connector(syn), i))
    ]

    eqs = cat(oldeqs, neweqs, dims = 1)
    systems = cat(ModelingToolkit.get_systems(group), synapse, dims = 1)

    return connected = ODESystem(eqs, t, [], [];
        name = name,
        systems = systems)
end

