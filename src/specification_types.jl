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

- voltage equation 
"""
struct BasicNeuron{G <: Geometry, D <: Dict, C <: FlowChannel} <: Neuron{T}
    geometry::G
    defaults::D
    channels::Vector{C}
    hooks::Int64
    name::Symbol
    function BasicNeuron{G,D,C}(g::G, d::D, c::Vector{C}, h::Int64, n::Symbol)
        T = Iterators.flatten(flows.(channels)) |> unique
        new{G,D,C,T}(g, d, c, h, n)
    end
end

# function BasicNeuron(g,d,c,h,n) 
#     T = Iterators.flatten(flows.(channels)) |> unique
#     return BasicNeuron(g,d,c,h,n,T)
# end

abstract type Geometry end 
abstract type FlowChannel end 
abstract type Neuron end



function StandardVoltageUpdate(::LiuSoma, neuron, flux)
    
end
