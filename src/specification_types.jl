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
struct BasicNeuron{G<:Geometry,D<:Dict,C<:FlowChannel} <: Neuron
    geometry::G
    defaults::D
    channels::Vector{C}
    hooks::Int64
    name::Symbol
end


LiuNeuron(d, vc, h, name) = BasicNeuron(NoGeometry(), d, vc, h, name)


function (b::BasicNeuron)()
    
    
    
    tracked = vcat(Iterators.flatten(sensedvars.(b.channels)) |> unique)
    tracked_vars = instantiate_variables(tracked)
    
    syns = [@variables $el(t) for el in [Symbol(:Isyn, i) for i = 1:hooks]]
    my_sum(syns) = hooks == 0 ? sum(Num.(syns)) : sum(reduce(vcat, syns))
    chs = [ch() for ch in b.channels]


    D(V)
    # tracked_fluxes = map(tracked) do thing
    #     if thing == :V
    #         sum(chs) do ch
    #             ModelingToolkit.getvar(ch.sys, currents )
    #         end        
    # end

end

function StandardVoltageDynamics(n::Neuron, flux)
    C = get_capacitance(n.geometry)
end

function get_capacitance(::NoGeometry)

end


abstract type FlowDynamics{F} end

struct BasicVoltageDynamics <: FlowDynamics{Voltage} end
"""
currents should include synaptic current
"""
function (b::BasicVoltageDynamics)(n::Neuron, V, currents, input)
    C = get_capacitance(n.geometry)
    return (1 / C) * (currents + input)
end
struct LiuCalciumDynamics{T<:Number} <: FlowDynamics{Calcium}
    τCa::T
    Ca∞::T
end

function (l::LiuCalciumDynamics)(n::Neuron, Ca, currents)
    C = get_capacitance(n.geometry)
    (1 / l.τCa) * (-Ca + Ca∞ + (f * area * summed_calcium_flux / C))
end

# defaults = Dict()
# neur = BasicNeuron(NoGeometry(20.), defaults, channels, 0, :AB)

