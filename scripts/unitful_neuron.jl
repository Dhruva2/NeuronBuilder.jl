using ModelingToolkit, OrdinaryDiffEq, Plots, NeuronBuilder, Unitful
import Unitful: mV, mS, cm, mm, nA, mA, µA, ms, nF, μM

const Cₘ = Cm() # specific capacitance cₘ is a biological constant (around) 10 nF/mm^2

Liu_conv = Cₘ

channels = [UnitLiu.Na(700.0 * Liu_conv), UnitLiu.CaS(4.0 * Liu_conv), UnitLiu.CaT(2.0 * Liu_conv), UnitLiu.Ka(50.0 * Liu_conv), UnitLiu.KCa(40.0 * Liu_conv),
    UnitLiu.Kdr(70.0 * Liu_conv), UnitLiu.H(0.03 * Liu_conv)]

#params = Params(Dict(:τCa => 200.0ms, :Ca∞ => 0.05µM, :Iapp => 0.0nA / mm^2))

#defaults doesn't have units yet
τCa = 200.0
Ca∞ = 0.05
fxarea = 14.96 * 0.0628

defaults = Dict(Voltage => BasicVoltageDynamics(),
    Calcium => Prinz.CalciumDynamics(τCa, Ca∞, fxarea),
    Reversal{Calcium} => Prinz.CaReversalDynamics())

ics = ICS(Dict(Voltage => -60.0mV, Calcium => 0.05μM))

reversals = Reversals(Dict(Reversal{Sodium} => 50.0mV,
    Reversal{Potassium} => -80.0mV,
    Reversal{Leak} => -50.0mV,
    Reversal{Proton} => -20.0mV,
    Reversal{Calcium} => 0.0mV))

somatic_parameters = merge(ics, reversals)

#somatic_parameters has to be converted back to subtype of real OR change BasicNeuron struct to accept Quantities 
b = BasicNeuron(NoGeometry(Cm), defaults, somatic_parameters, channels, :test_UnitLiu)

neur = b()

prob = ODEProblem(neur, [], (0.0, 5000.0), [])

sol = solve(prob, AutoTsit5(Rosenbrock23()))
plot(sol)