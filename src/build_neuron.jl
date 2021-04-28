
function build_channel(ch::IonChannel; name = get_name(ch))
    @parameters t
    D = Differential(t)
    @variables V(t) Ca(t)
    eqs, states, parameters, current, defaultmap = channel_dynamics(ch, V, Ca, D, t)
    return ODESystem(eqs, t, [V, Ca, states...], [parameters...]; 
                                            observed=current, 
                                            name = name, 
                                            defaults = defaultmap)
end


function build_channel(syn::Synapse; name = get_name(syn))
    @parameters t
    D = Differential(t)
    @variables Vpre(t) Vpost(t)
    eqs, states, parameters, current, defaultmap = channel_dynamics(syn, Vpre, Vpost, D, t)
    return ODESystem(eqs, t, [Vpre, Vpost, states...], [parameters...];    observed = current,
                        name = name,
                        defaults = defaultmap)
end


# 10 nF/mm² for default specific membrane capacitance
function build_neuron(channels;  capacitance = 10, V_init = -60, Ca_init = 0.05, name = :unidentified_neuron, num_inputs = 0)

    @parameters t cm
    states = @variables V(t) Ca(t)
    syns = [Num(Variable{Symbolics.FnType{Tuple{Any},Real}}(Symbol(:Isyn, i)))(t) for i in 1:num_inputs]

    D = Differential(t)
    channel_systems = build_channel.(channels)
    
    summed_membrane_currents = sum(ionic_current(ch, ch_sys) for (ch, ch_sys) in zip(channels, channel_systems))

    summed_calcium_flux = sum(calcium_current(ch, ch_sys) for (ch, ch_sys) in zip(channels, channel_systems))

    connections = cat(
                    [voltage_hook(V, ch, chs) for (ch, chs) in zip(channels, channel_systems)],
                    [calcium_hook(Ca, ch, chs) for (ch, chs) in zip(channels, channel_systems)], 
                    dims=1
                    )  

    diffeqs =   cat(
                    [D(V) ~ (summed_membrane_currents + sum(syns))/cm],
                [D(Ca) ~ (summed_calcium_flux*.0939 - Ca + .05)/200],
                dims=1)               

    eqs = cat(connections, diffeqs; dims=1)
    eqs = filter(x -> !(x == (Num(0) ~ Num(0))), eqs)

    all_states = [states..., syns...]

    neur = ODESystem(eqs, t; 
                    systems = channel_systems,
                    defaults = [V => V_init, Ca => Ca_init, cm => capacitance],
                    name = name
                    )

    return neur
end


function build_group(arr_of_neurons; name = :group)
    @parameters t 
    eqs = Array{Equation,1}(undef,0)
    ODESystem(eqs, t, [], []; name=name, systems = arr_of_neurons)
end


"""
Issue: need to do all connections in one ḡChol

"""
function add_connection(group, pre, post, syn::Synapse; name = ModelingToolkit.getname(group), i=1)

    @parameters t
    prename, postname = ModelingToolkit.getname.([pre, post])
    synapse = build_channel(syn; name = Symbol(prename, :to, postname, get_name(syn)))

    oldeqs = ModelingToolkit.get_eqs(group)
    neweqs = [
            get_pre(syn, pre) ~ my_pre(syn, synapse),
            get_post(syn, post) ~ my_post(syn, synapse),
            ionic_current(syn, synapse) ~ getproperty(post, Symbol(post_connector(syn), i))
    ]
    
    eqs = cat(oldeqs, neweqs, dims=1)
    systems = cat(ModelingToolkit.get_systems(group), synapse, dims=1)

    return connected = ODESystem(eqs, t, [], []; 
                                name=name, 
                                systems = systems)
end

