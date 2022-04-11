abstract type Specification end

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


# Iterators.flatten(reversals.(channels)) |> unique


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
struct BasicNeuron{G<:Geometry,C<:FlowChannel} <: Neuron
    geometry::G
    dynamics::Dict{DataType,SpeciesDynamics}
    somatic_parameters::Dict{DataType,Float64}
    channels::Vector{C}
    name::Symbol
end


LiuNeuron(d, vc, h, name) = BasicNeuron(NoGeometry(), d, vc, h, name)


function (b::BasicNeuron)(hooks::Integer)

    # if F ∉ b.dynamics, return x->0. else return Fdynamics(neur, var, rhs)
    dynamics(F, neur, vars, varnames, RHS) =
        get(neur.dynamics, F) do
            f(w, x, y, z) = Num(0.0) ~ Num(0.0)
        end(neur, vars, varnames, RHS)


    # track every species sensed by the connected channels
    tracked = vcat(Voltage, b.channels .|> sensed |> Iterators.flatten |> unique)

    # build state variables for each of these tracked species 
    tracked_vars = tracked .|> shorthand_name |> instantiate_variables


    syns = [@variables $el(t) for el in [Symbol(:Isyn, i) for i = 1:hooks]]
    my_sum(syns) = hooks == 0 ? sum(Num.(syns)) : sum(reduce(vcat, syns))
    chs = [ch() for ch in b.channels]


    tracked_fluxes = map(tracked) do thing
        sum(chs) do ch
            get_actuator(ch, thing)
        end
    end

    ## UGLY, make better. maybe add generality for syns:. Like hooks = {Voltage, 3}
    tracked_fluxes[findfirst(x -> x == Voltage, tracked)] += my_sum(syns)

    ## ie hook connections between channel sensors and tracked variables. All LHS are variables. Some RHS are parameters.
    outward_connections = reduce(vcat, map(tracked, tracked_vars) do species, variable
        [variable ~ get_sensor(ch, species) for ch in chs if (get_sensor(ch, species) !== nothing)]
    end)

    outward_species = reduce(vcat, map(tracked) do species
        [species for ch in chs if (get_sensor(ch, species) !== nothing)]
    end
    )

    """
    1.  Get rid of outward connections where RHS is a parameter (i.e. static reversals)
    2.  Save these indices 
    3.  For each of the parameters on the RHS of these indices:
        a. Remove the corresponding variable from tracked variables
        b. Make a dictionary with the parameter going to its default

    For ECa: 
    - delete the algebraic equation in the cas and cat channels
    - Add a SpeciesDynamics
"""

    ## ie add currents and inputs to the dynamics of tracked variables
    inward_connections = [
        dynamics(species, b, tracked_vars, tracked, flux) for (var, flux, species) in zip(tracked_vars, tracked_fluxes, tracked)
    ] |> y -> filter(x -> !isnothing(x), y)

    incorrect_indices = findall(outward_connections) do el
        ModelingToolkit.isparameter(el.rhs)
    end

    outward_connection_indices = setdiff(1:length(outward_connections), incorrect_indices)

    parameter_defaults = Dict(outward_connections[el].rhs => b.somatic_parameters[outward_species[el]] for el in unique(incorrect_indices))

    param_indices = findall(x -> x ∈ outward_species[incorrect_indices], tracked)
    state_indices = setdiff(1:length(tracked), param_indices)

    state_vars = tracked_vars[state_indices]
    state_defaults = Dict(tracked_vars[el] => b.somatic_parameters[tracked[el]] for el in state_indices)

    return ComponentSystem(b, ODESystem(
        vcat(inward_connections, outward_connections[outward_connection_indices]),
        t;
        systems=[ch.sys for ch in chs],
        defaults=merge(state_defaults, parameter_defaults),
        name=b.name
    ) |> structural_simplify)
end

# |> ch -> filter(x -> x.rhs !== nothing, ch)
# StandardVoltageDynamics(b, (tracked_vars) .~ tracked_fluxes)


# function StandardVoltageDynamics(n::Neuron, V, flux)
#     C = get_capacitance(n.geometry)
# end





# defaults = Dict()
# neur = BasicNeuron(NoGeometry(20.), defaults, channels, 0, :AB)

"""
final decisions:
tracked variables includes all the reversals

Neurons have dynamics for all the reversals. These default to zero for the BasicNeuronReversal.

At the end of system construction, take every state that has zero dynamics and turn it into a parameter.

Inputs are not defined. If you want inputs, make a channel. 
"""


"""
channel defs:
if reversal is a parameter, then hook it to something 
if reversal is already taken (either by being a variable with an equation, or having a default) then 

external parameters are:
reversals (where they are parameters)

currently for new setup, ECa(t) is defined in the channel
    more generally it will depend upon neural parameters
    so really it is a sensor

Really: channel should ask neuron what the type of the reversal is, and then hook param ̃ param or variable ˜ variable 

Or there should be an automatic hooking of reversals that works for params or vars 


"""

"""
Dealing with synapses:

it's just a layer on top of the same type.

Right now: in neuron:
build each ion channel()
soma tracks union of sensed variables from each ion channel
tracked dynamics use each channel for calculation 

NeuronGroup type. 
in neuron group
build each neuron, with functions for nodes and edges 

"""