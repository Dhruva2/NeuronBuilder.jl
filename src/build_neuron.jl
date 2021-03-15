
function build_channel(ch::IonChannel)
    @parameters t
    D = Differential(t)
    @variables V(t) Ca(t)
    eqs, states, parameters, current, u0map, pmap = channel_dynamics(ch, V, Ca, D, t)
    return ODESystem(eqs, t, [V, Ca, states...], [parameters...]; 
                                            observed=current, 
                                            name = get_name(ch), 
                                            default_u0 = u0map, 
                                            default_p = pmap)    
end


function build_channel(syn::Synapse)
    @parameters t
    D = Differential(t)
    @variables Vpre(t) Vpost(t)
    eqs, states, parameters, current, u0map, pmap = channel_dynamics(syn, Vpre, Vpost, D, t)
    return ODESystem(eqs, t, [Vpre, Vpost, states...], [parameters...];    observed = current,
                        name = get_name(syn),
                        default_u0 = u0map,
                        default_p = pmap)
end



function build_neuron(channels; synapses = [], V_init = -60, Ca_init = 0.05)

    @parameters t 
    states = @variables V(t) Ca(t)
    D = Differential(t)
    channel_systems = build_channel.(channels)
    synapse_systems = build_synapse.(synapses)

    summed_membrane_currents = sum(ionic_current(ch, ch_sys) for (ch, ch_sys) in zip(channels, channel_systems))

    if length(synapases) > 0
        summed_synaptic_currents = sum(synaptic_current(syn, syn_sys) for (syn, syn_sys) in zip(synapses, synapse_systems))
    else
        summed_synaptic_currents = Num(0.)
    end

    summed_calcium_flux = sum(calcium_current(ch, ch_sys) for (ch, ch_sys) in zip(channels, channel_systems))

    alg_connections = cat(
                    [V ~ el.V for el in channel_systems],
                    [Ca ~ el.Ca for el in channel_systems],
                    [V ~ el.Vpost for el in synapse_systems],
                    dims=1
                    )  

    diffeqs =   cat(
                [D(V) ~ summed_membrane_currents + summed_synaptic_currents],
                [D(Ca) ~ summed_calcium_flux],
                dims=1)               

    eqs = cat(alg_connections, diffeqs; dims=1)
    all_states = [states...]

    neur = ODESystem(eqs, t, all_states, []; 
                    systems = channel_systems,
                    default_u0 = [V => V_init, Ca => Ca_init]
                    )
end


function ff(x)
    if x == 1
        return nothing
    else
        return x
end