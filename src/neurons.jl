
# struct Neuron{C<:Component, S<:Soma,OS<:ODESystem}
#     soma::S
#     channels::Vector{C}
#     hooks::Int64
#     name::Symbol
#     ODESystem::OS
# end

# Neuron(a, b, c, d) = Neuron(a, b, c, d, build_neuron(a, b, c, d))
# Neuron(s, t) = Neuron(s, t, 0, :unidentified_neuron)
# Neuron(n::Neuron, new_hooks::Integer) = Neuron(n.soma, n.channels, new_hooks, n.name)


struct BasicVoltageDynamics <: SpeciesDynamics{Voltage} end

function (b::BasicVoltageDynamics)(n::Neuron, vars, varnames, flux)
    Cₘ, = get_parameters(b)
    V = vars[findfirst(x -> x == Voltage, varnames)]
    return D(V) ~ (1 / Cₘ) * (flux)
end
get_parameters(::BasicVoltageDynamics) = @parameters Cₘ
default_params(v::BasicVoltageDynamics, n::Neuron, vars, varnames) = Dict(get_parameters(v)... => capacitance(n.geometry))



"""
building neuron 
- geometric parameters are unique and should be in compartment. also optional (or = 1)
- same for electrical parameters ie capacitance 

- instead of number of hooks, why not connection as a type of channel?

- reversals are defined by the ion channels to which it is connected. so not in type signature
- same for currents
- update equations: 
- voltage equation 
"""

struct EmptyNeuron{F <: Number} <: Neuron
    somatic_parameters::Dict{DataType, F}
end

EmptyNeuron() = EmptyNeuron(
    Dict(Voltage => -60.)
)
dynamics(::EmptyNeuron) = Dict{DataType, SpeciesDynamics}()


struct BasicNeuron{G<:Geometry,C<:FlowChannel, F <: Real} <: Neuron
    geometry::G
    dynamics::Dict{DataType,SpeciesDynamics}
    somatic_parameters::Dict{DataType,F}
    channels::Vector{C}
    name::Symbol
end


LiuNeuron(d, vc, h, name) = BasicNeuron(NoGeometry(), d, vc, h, name)



function (b::BasicNeuron)(hooks::Integer)
    # DO LATER: to make it easier for the user, add an extra input which is var (the var in question). so they dont have to find eg voltage by searching through vars and varnames
    # add b as input to channels so they can sense whether their inputs are parmeters or not, using b.dynamics

    # e.g. species = Voltage or species = Potassium
    has_dynamics(species) = haskey(b.dynamics, species)

    # track union of things sensed by the connected channels
    tracked_names = vcat(Voltage, b.channels .|> sensed |> Iterators.flatten |> unique)

    state_indices = findall(has_dynamics, tracked_names)
    param_indices = findall(!has_dynamics, tracked_names)

    # build state/param ModelingToolkit variables for each of these tracked species 
    tracked = zeros(Num, length(tracked_names))
    
    tracked[state_indices] = reduce(vcat,
        tracked_names[state_indices] .|> shorthand_name
    ) |> instantiate_variables

    tracked[param_indices] = reduce(vcat,
        tracked_names[param_indices] .|> shorthand_name
    ) |> instantiate_parameters


    syns = [@variables $el(t) for el in [Symbol(:Isyn, i) for i = 1:hooks]]
    my_sum(syns) = hooks == 0 ? sum(Num.(syns)) : sum(reduce(vcat, syns))
    !(hooks == 0) && (syns = reduce(vcat, syns))
    chs = [ch(b) for ch in b.channels]


    tracked_fluxes = map(tracked_names[state_indices]) do thing
        sum(chs) do ch
            get_actuator(ch, thing)
        end
    end

    ## UGLY, make better. maybe add generality for syns:. Like hooks = {Voltage, 3}
    tracked_fluxes[findfirst(x -> x == Voltage, tracked_names)] += my_sum(syns)

    ## ie hook connections between channel sensors and tracked variables. 
    outward_connections = reduce(vcat, map(tracked_names, tracked) do species, variable
        [variable ~ get_sensor(ch, species) for ch in chs if (get_sensor(ch, species) !== nothing)]
    end)

    outward_species = reduce(vcat, map(tracked_names) do species
        [species for ch in chs if (get_sensor(ch, species) !== nothing)]
    end
    )

    outward_connection_indices = findall(outward_connections) do el
        !ModelingToolkit.isparameter(el.rhs)
    end

    ## ie add currents and inputs to the dynamics of tracked variables
    inward_connections = [
        b.dynamics[species](b, tracked, tracked_names, flux) for (flux, species) in zip(tracked_fluxes, tracked_names[state_indices])
    ]

    # default initial conditions for neural state variables
    state_defaults = Dict(tracked[el] => b.somatic_parameters[tracked_names[el]] for el in state_indices)

    # default values for hooked parameter values in channels (their internal states/parameters already have defaults)
    outward_param_indices = findall(outward_connections) do el
        ModelingToolkit.isparameter(el.rhs)
    end
    parameter_defaults = Dict(outward_connections[el].rhs => b.somatic_parameters[outward_species[el]] for el in unique(outward_param_indices))

    # if there are any internal state variables from the b.dynamics laws, store and default them.
    somatic_states = get_from(b.dynamics, get_states)
    somatic_state_defaults = mapreduce(merge, values(b.dynamics)) do s
        default_states(s, b, tracked, tracked_names)
    end

    # if there are any internal parameter variables from the b.dynamics laws, store and default them.
    somatic_params = get_from(b.dynamics, get_parameters)
    somatic_param_defaults = mapreduce(merge, values(b.dynamics)) do s
        default_params(s, b, tracked, tracked_names)
    end

    sys = ODESystem(
        vcat(inward_connections, outward_connections[outward_connection_indices]),
        t,
        vcat(tracked[state_indices], somatic_states, syns),
        somatic_params
        ;
        systems=[ch.sys for ch in chs],
        defaults=merge(state_defaults, parameter_defaults, somatic_state_defaults, somatic_param_defaults),
        name=b.name
    )
    return ComponentSystem(b, (hooks == 0) ? structural_simplify(sys) : sys)
end

