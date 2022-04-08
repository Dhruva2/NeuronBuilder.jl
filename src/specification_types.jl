abstract type Specification end

struct Neuron{S<:Soma,OS<:ODESystem}
    soma::S
    channels::Vector{Component}
    hooks::Int64
    name::Symbol
    ODESystem::OS
end


Neuron(a, b, c, d) = Neuron(a, b, c, d, build_neuron(a, b, c, d))
Neuron(s, t) = Neuron(s, t, 0, :unidentified_neuron)
Neuron(n::Neuron, new_hooks::Integer) = Neuron(n.soma, n.channels, new_hooks, n.name)


