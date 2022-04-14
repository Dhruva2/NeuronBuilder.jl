"""
Prinz bursting neuron script
"""

using NeuronBuilder, ModelingToolkit, OrdinaryDiffEq, Plots

const Area = 0.0628 # Prinz/Liu 0.0628 mm2
const Cm = 10.0 # specific capacitance cₘ is a biological constant (around) 10 nF/mm^2

Prinz_conv = Cm / Area

channels = [Prinz.Na(100 * Prinz_conv), Prinz.CaS(6 * Prinz_conv), Prinz.CaT(2.5 * Prinz_conv), Prinz.H(0.01 * Prinz_conv),
    Prinz.Ka(50 * Prinz_conv), Prinz.KCa(5 * Prinz_conv), Prinz.Kdr(100 * Prinz_conv), Prinz.leak(0.0)]

τCa = 200.0
Ca∞ = 0.05

defaults = Dict(Voltage => BasicVoltageDynamics(),
    Calcium => Prinz.CalciumDynamics(τCa, Ca∞, 14.96 * 0.0628),
    Reversal{Calcium} => Prinz.CaReversalDynamics())

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