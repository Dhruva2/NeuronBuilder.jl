
function build_channel(ca::T) where T<:Channel
    @parameters t
    D = Differential(t)
    @variables V(t) Ca(t)
    eqs, states, parameters, current, u0map, pmap = channel_dynamics(ca, V, Ca, D, t)
    return ODESystem(eqs, t, [V, Ca, states...], [parameters...]; 
                                            observed=current, 
                                            name = get_name(ca), 
                                            default_u0 = u0map, 
                                            default_p = pmap)    
end




function build_neuron(channels; number_of_inputs=0)

    @parameters t 
    states = @variables V(t) Ca(t)
    @variables INa(t) ICaS(t) ICaS_Ca(t)
    D = Differential(t)
    all_systems = build_channel.(channels)
   
    alg_connections = cat(
                    [V ~ el.V for el in all_systems],
                    [Ca ~ el.Ca for el in all_systems],
                    dims=1
                    )  

    diffeqs =   cat(
                [D(V) ~ sum(voltage_current(ch, ch_sys) for (ch, ch_sys) in zip(channels, all_systems))],
                [D(Ca) ~ sum(calcium_current(ch, ch_sys) for (ch, ch_sys) in zip(channels, all_systems))],
                dims=1)               

    eqs = cat(alg_connections, diffeqs; dims=1)
    all_states = [states...]

    neur = ODESystem(eqs, t, all_states, []; systems = all_systems)
end