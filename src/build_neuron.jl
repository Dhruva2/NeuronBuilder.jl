function (ch::IonChannel)()
    @variables V(t) Ca(t)
    eqs, states, parameters, current, defaultmap = channel_dynamics(ch, V, Ca)
    sys = ODESystem(eqs, t, [V, Ca, states...], [parameters...];
        observed = current, defaults = defaultmap, name = get_name(ch))
    return ComponentSystem(ch,sys)
end

function (ch::IonChannel)(reg::Bool, Ca_tgt, τg)
    if isLeak(ch) || (reg == false)
        ch()
    elseif (reg == true)
        @variables V(t) Ca(t)
        eqs, states, parameters, current, defaultmap = channel_dynamics(ch, V, Ca)
        sys = ODESystem(eqs, t, [V, states...], [parameters...];
            observed = current, defaults = defaultmap, name = get_name(ch))
        @parameters Ca_tgt, τg
        eqs, states, parameters, defaultmap = regul_dyn(ch, Ca, Ca_tgt, τg)
        @named regul_sys = ODESystem(eqs, t, [Ca, states...], parameters, defaults = defaultmap)
        sys = extend(sys, regul_sys)
    
        return ComponentSystem(ch, sys)
    end
end

function build_channel(syn::Synapse; name = get_name(syn))
    @variables Vpre(t) Vpost(t)
    eqs, states, parameters, current, defaultmap = channel_dynamics(syn, Vpre, Vpost)
    sys = ODESystem(eqs, t, [Vpre, Vpost, states...], [parameters...]; observed = current,
        name = name, defaults = defaultmap)
    return ComponentSystem(syn,sys)
end


function build_neuron(comp, channels; reg::Bool = false, name = :unidentified_neuron)

    neuron_parameters = @parameters Cₘ τCa Ca∞ Ca_tgt τg Iapp
    neuron_states = @variables V(t) Ca(t)
    syns = [@variables $el(t) for el in [Symbol(:Isyn, i) for i = 1:comp.hooks]]
    my_sum(syns) = comp.hooks == 0 ? sum(Num.(syns)) : sum(reduce(vcat, syns))

    channel_systems = [ch(reg, Ca_tgt, τg) for ch in channels]

    summed_membrane_currents = sum(ionic_current(ch, cs.sys) for (ch, cs) in zip(channels, channel_systems))

    summed_calcium_flux = sum(calcium_current(ch, cs.sys) for (ch, cs) in zip(channels, channel_systems))

    connections = cat(
        [voltage_hook(V, cs) for cs in channel_systems],
        [calcium_hook(Ca, cs) for cs in channel_systems],
        dims = 1
    )
    f = 0.94
    diffeqs = cat(
        [D(V) ~ (summed_membrane_currents + my_sum(syns) + Iapp) / Cₘ],
        [D(Ca) ~ (1 / τCa) * (-Ca + Ca∞ + (f * summed_calcium_flux))], #ionic currents are already carrying negative sign
        dims = 1
    )

    eqs = cat(connections, diffeqs; dims = 1)
    eqs = filter(x -> !(x == (Num(0) ~ Num(0))), eqs)

    channel_defs = Dict{Num,Float64}()
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
    synapse_sys = build_channel(syn; name = Symbol(prename, :to, postname, get_name(syn)))

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

