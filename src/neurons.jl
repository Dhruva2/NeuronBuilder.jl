dynamics(n::Neuron) = n.dynamics
has_dynamics(n::Neuron, species) = haskey(dynamics(n), species)

# defined for reversal and conductance since they belong to a neuron. not conductance (belongs to Link)
has_dynamics(n::Neuron, ::Type{Current{I}}) where {I} = true
# has_dynamics(n::Neuron, ::Type{Reversal{I}}) where {I} = haskey(dynamics(n), I)


struct ResetDynamics{T<:Number} <: SpeciesDynamics{Voltage}
    V_threshold::T
    V_reset::T
    function ResetDynamics(x::T, y::T) where {T<:Number}
        x < y ? error("Threshold lower value than reset.") : new{T}(x, y)
    end
end

get_parameters(::ResetDynamics) = @parameters Cₘ V_threshold V_reset
default_params(l::ResetDynamics, n::Neuron, vars, varnames) = Dict(
    get_parameters(l) .=> (capacitance(n.geometry), l.V_threshold, l.V_reset)
)

function (b::ResetDynamics)(n::Neuron, vars, varnames, flux)
    Cₘ, V_threshold, V_reset = get_parameters(b)
    V = vars[findfirst(x -> x == Voltage, varnames)]
    return D(V) ~ (1 / Cₘ) * (flux) #non-standard convention, sum of fluxes already has negative sign because of the (E-V) in currents
end

kwargs(::SpeciesDynamics, vars, varnames) = NamedTuple()
function kwargs(b::ResetDynamics, vars, varnames) 
    V = vars[findfirst(x -> x == Voltage, varnames)]
    reset = [V ~ b.V_threshold] => [V ~ b.V_reset]
    return (:continuous_events => reset,)
end
struct EmptyNeuron{F<:Number} <: Neuron
    somatic_parameters::Dict{DataType,F}
end

EmptyNeuron() = EmptyNeuron(
    Dict(Voltage => -60.0)
)
dynamics(::EmptyNeuron) = Dict{DataType,SpeciesDynamics}()


struct BasicNeuron{G<:Geometry,C<:FlowChannel,F<:Real,S<:SpeciesDynamics} <: Neuron
    geometry::G
    dynamics::Dict{DataType,S}
    somatic_parameters::Dict{DataType,F}
    channels::Vector{C}
    name::Symbol
end


"""
building neuron 
- geometric and electrical (ie capacitance) parameters are unique and should be defined by user
- number of hooks has to be pre-specified
- reversals are defined by the ion channels to which it is connected, same for currents
- channels which have a PlasticityRule (need extra sensors) are treated inside b the same way as those that don't 
- update equations: voltage equation gets added synapses...for now. TODO Make synapses like channels

# DO LATER: to make it easier for the user, add an extra input which is var (the var in question). so they dont have to find eg voltage by searching through vars and varnames
# add b as input to channels so they can sense whether their inputs are parmeters or not, using b.dynamics
"""

function (b::BasicNeuron)(; incoming_connections::Union{Integer, Bool} = false)
    #shared -> hooks
    # e.g. species = Voltage or species = Potassium
    has_dynamics(species) = haskey(b.dynamics, species)
    # track union of things sensed by the connected channels

    tracked_names = vcat(Voltage, b.channels .|> sensed |> el -> filter(!isnothing, el) |> Iterators.flatten |> unique)
    tracked_names = vcat(
        Voltage,
        b.channels .|> sensed |> Iterators.flatten |> unique,
        [keys(b.dynamics)...]
    ) |> unique!

    state_indices = findall(has_dynamics, tracked_names)
    param_indices = findall(!has_dynamics, tracked_names)

    # build state/param ModelingToolkit variables for each of these tracked species 
    tracked = zeros(Num, length(tracked_names))

    tracked[state_indices] .= reduce(vcat,
        tracked_names[state_indices] .|> shorthand_name
        |> instantiate_variables)

    tracked[param_indices] .= reduce(vcat,
        tracked_names[param_indices] .|> shorthand_name
        |> instantiate_parameters)

    syns = [@variables $el(t) for el in [Symbol(:Isyn, i) for i = 1:incoming_connections]]
    my_sum(syns) = incoming_connections == 0 ? 0.0 : sum(reduce(vcat, syns))
    !(incoming_connections == 0) && (syns = reduce(vcat, syns))
    chs = [ch(b) for ch in b.channels]

    tracked_fluxes = map(tracked_names[state_indices]) do thing
        sum(zip(b.channels, chs)) do (channel, channel_sys)
            get_actuator(channel, channel_sys, thing)
        end |> Num
    end 

    # UGLY, make better. maybe add generality for syns:. Like hooks = {Voltage, 3}
    tracked_fluxes[findfirst(x -> x == Voltage, tracked_names)] += my_sum(syns)

    # hook connections between channel sensors and tracked variables. 
    outward_connections = reduce(vcat, map(tracked_names, tracked) do species, variable
        [(variable ~ get_sensor(ch, ch_sys, species))::Equation for (ch, ch_sys) in zip(b.channels, chs) if (get_sensor(ch, ch_sys, species) !== nothing)]::Vector{Equation}
    end)

    outward_species = reduce(vcat, map(tracked_names) do species
        [species for (ch, ch_sys) in zip(b.channels, chs) if (get_sensor(ch, ch_sys, species) !== nothing)]
    end
    )

    outward_connection_indices = findall(outward_connections) do el
        !ModelingToolkit.isparameter(el.rhs)
    end

    # add currents and inputs to the dynamics of tracked variables
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
        systems=chs,
        defaults=merge(state_defaults, parameter_defaults, 
        somatic_state_defaults, somatic_param_defaults),
        ((kwargs(el, tracked, tracked_names) for el in values(b.dynamics))...)...,
        name=b.name,
    )
    if typeof(incoming_connections) == Bool && !incoming_connections
        return structural_simplify(sys)
    elseif typeof(incoming_connections) <: Integer
        return sys
    end


end

export get_from, get_states, default_states, default_params, get_parameters
