"""
Liu bursting neuron script
"""

using NeuronBuilder, ModelingToolkit, OrdinaryDiffEq, Plots

const Area = 0.0628 # Prinz/Liu 0.0628 mm2
const Cm = 10.0 # specific capacitance cₘ is a biological constant (around) 10 nF/mm^2

Liu_conv = Cm

channels = [Liu.Na(700.0 * Liu_conv), Liu.CaS(4.0 * Liu_conv), Liu.CaT(2.0 * Liu_conv), Liu.Ka(50.0 * Liu_conv), Liu.KCa(40.0 * Liu_conv),
    Liu.Kdr(70.0 * Liu_conv), Liu.H(0.03 * Liu_conv), Liu.leak(0.01 * Liu_conv)]

τCa = 20.0
Ca∞ = 0.05

defaults = Dict(Voltage => BasicVoltageDynamics(),
    Calcium => LiuCalciumDynamics(τCa, Ca∞, 14.96 * 0.0628), 
    Reversal{Calcium} => LiuCaReversalDynamics())

somatic_parameters = Dict(
    Reversal{Sodium} => 50.0,
    Reversal{Potassium} => -80.0,
    Reversal{Leak} => -50.0,
    Reversal{Proton} => -20.0,
    Voltage => -60.0,
    Calcium => 0.05,
    Reversal{Calcium} => 0.0)

b = BasicNeuron(NoGeometry(Cm), defaults, somatic_parameters, channels, :test_Liu)

neur = b(0)

prob = ODEProblem(neur.sys, [], (0.0, 5000.0), [])

sol = solve(prob, Tsit5())
plot(sol)
