"""
LIF neuron script
"""

using NeuronBuilder, ModelingToolkit, OrdinaryDiffEq, Plots
struct leak{D} <: FlowChannel(Leak)
    g::D
end

function (ch::leak)(n::Neuron; name=get_name(ch))
    I, E, V = instantiate_hooks(n, ch, currents, reversals, voltage)
    g, = instantiate_parameters(ch, conductances)

    eqs = [I ~ g * (E - V)]

    states, params = vardivide(V, E, I, g)
    defaultmap = [g => ch.g]

    return ODESystem(eqs, t, states, params; defaults=defaultmap, name=name)
end
struct Input{F} <: FlowChannel{Tuple{Nothing}, Tuple{Leak}}
    f::F
end

function (ch::Input)(neur::Neuron; name=get_name(ch))
    I, V = instantiate_hooks(neur, ch, currents, voltage)

    eqs = [I ~ ch.f(t)]
    states, params = vardivide(V, I)
    return ODESystem(eqs, t, states, params; defaults=[I => ch.f(0)], name)
end

input_current = Input(t -> 10.)
channels = [leak(10.), input_current]

defaults = Dict{DataType, SpeciesDynamics}(Voltage => ResetDynamics(1., 0.))

Cm = 10.0
somatic_parameters = Dict(
    Voltage => 0.0,
    Reversal{Leak} => 0.0
    )

b = BasicNeuron(NoGeometry(Cm), defaults, somatic_parameters, channels, :test_LIF)

neur = b()

prob = ODEProblem(neur, [], (0.0, 50.0), [])

sol = solve(prob, Tsit5(); abstol = 1e-10, reltol = 1e-10)
plot(sol)
