"""
Observing a HH neuron and its parameters
"""

using NeuronBuilder, ModelingToolkit, OrdinaryDiffEq, Plots

Cₘ = 1.0
channels = [HH.Na(120),HH.K(36),HH.leak(0.3)]
dynamics = Dict(
    Voltage => BasicVoltageDynamics())
somatic_parameters = Dict(
    Reversal{Sodium} => 55.0,
    Reversal{Potassium} => -77.0,
    Reversal{Leak} => -54.4,
    Voltage => -60.0)

hh = BasicNeuron(NoGeometry(Cₘ), dynamics, somatic_parameters, channels, :HH) #point neuron
neur = hh()
