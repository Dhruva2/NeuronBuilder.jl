"""
Prinz bursting neuron script
"""

using NeuronBuilder, ModelingToolkit, OrdinaryDiffEq, Plots

const Area = 0.0628 # Prinz/Liu 0.0628 mm2
const Cm = 10.0 # specific capacitance cₘ is a biological constant (around) 10 nF/mm^2

Prinz2_conv = Cm / Area

channels = [Prinz2.Na(100 * Prinz2_conv), Prinz2.CaS(6 * Prinz2_conv), Prinz2.CaT(2.5 * Prinz2_conv), Prinz2.H(0.01 * Prinz2_conv),
    Prinz2.Ka(50 * Prinz2_conv), Prinz2.KCa(5 * Prinz2_conv), Prinz2.Kdr(100 * Prinz2_conv), Prinz2.leak(0.0)]

τCa = 200.0
Ca∞ = 0.05

defaults = Dict(Voltage => BasicVoltageDynamics(),
    Calcium => PrinzCalciumDynamics(τCa, Ca∞, 14.96 * 0.0628),
    Reversal{Calcium} => PrinzCaReversalDynamics())

somatic_parameters = Dict(
    Reversal{Sodium} => 50.0,
    Reversal{Potassium} => -80.0,
    Reversal{Leak} => -50.0,
    Reversal{Proton} => -20.0,
    Voltage => -60.0,
    Calcium => 0.05,
    Reversal{Calcium} => 0.0)

b = BasicNeuron(NoGeometry(Cm), defaults, somatic_parameters, channels, :test_Prinz)

neur = b(0)

prob = ODEProblem(neur.sys, [], (0.0, 5000.0), [])

sol = solve(prob, AutoTsit5(Rosenbrock23()))
plot(sol)