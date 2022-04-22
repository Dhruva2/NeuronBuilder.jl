"""
LIF neuron script
"""

using NeuronBuilder, ModelingToolkit, OrdinaryDiffEq, Plots

struct leak{D<:Real} <: FlowChannel(Leak)
    gLeak::D
end

function (ch::leak)(n::Neuron; name=get_name(ch))
    I, E, V = instantiate_hooks(n, ch, currents, reversals, voltage)
    g, = instantiate_parameters(ch, conductances)

    eqs = [I ~ g * (E - V)]

    states, params = vardivide(V, E, I, g)
    defaultmap = [g => ch.gLeak]

    return ODESystem(eqs, t, states, params; defaults=defaultmap, name=name)
end

FlowChannel() = FlowChannel{Actuators <: Tuple}
struct Input <: FlowChannel{Tuple{Nothing}, Tuple{PseudoIon}} end

function (ch::Input)(n::Neuron; name = get_name(ch))
    I, V = instantiate_hooks(n, ch, currents, voltage)

    eqs = [I ~ 5.]

    states, params = vardivide(V, I)

    return ODESystem(eqs, t, states, params; name=name)
end

channels = [leak(10.), Input()]

defaults = Dict{DataType, SpeciesDynamics}(Voltage => ResetDynamics(1., 0.))

Cm = 10.0
somatic_parameters = Dict(
    Voltage => 0.0,
    Reversal{Leak} => 0.0
    )

b = BasicNeuron(NoGeometry(Cm), defaults, somatic_parameters, channels, :test_LIF)

neur = b()

prob = ODEProblem(neur, [], (0.0, 50.0), [])

sol = solve(prob, Tsit5())
plot(sol)
