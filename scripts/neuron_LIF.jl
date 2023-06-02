"""
LIF neuron script
"""

using ModelingToolkit, OrdinaryDiffEq, Plots


input_current = NeuronBuilder.HodgkinHuxley.Input(t -> 10.0)
channels = [NeuronBuilder.Liu.leak(10.), input_current]
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
