"""
Liu bursting neuron script
"""

using ModelingToolkit, OrdinaryDiffEq, Plots

const Cm = 10.0 # specific capacitance cₘ is a biological constant (around) 10 nF/mm^2
Liu_conv = Cm

channels = [Liu.Na(700.0 * Liu_conv), Liu.CaS(4.0 * Liu_conv), Liu.CaT(2.0 * Liu_conv), Liu.Ka(50.0 * Liu_conv), Liu.KCa(40.0 * Liu_conv),
    Liu.Kdr(70.0 * Liu_conv), Liu.H(0.03 * Liu_conv), Liu.leak(0.01 * Liu_conv)]

τCa = 20.0
Ca∞ = 0.05
fxarea = 14.96 * 0.0628

dynamics = Dict(Voltage => BasicVoltageDynamics(),
    Calcium => Liu.CalciumDynamics(τCa, Ca∞, fxarea),
    Reversal{Calcium} => Liu.CaReversalDynamics())

somatic_parameters = Dict(
    Reversal{Sodium} => 50.0,
    Reversal{Potassium} => -80.0,
    Reversal{Leak} => -50.0,
    Reversal{Proton} => -20.0,
    Voltage => -60.0,
    Calcium => 0.05,
    Reversal{Calcium} => 0.0)

b = BasicNeuron(NoGeometry(Cm), dynamics, somatic_parameters, channels, :test_Liu) #point neuron

neur = b()

prob = ODEProblem(neur, [], (0.0, 5000.0), []; jac=true)

sol = @time solve(prob, Tsit5())
plot(sol)
