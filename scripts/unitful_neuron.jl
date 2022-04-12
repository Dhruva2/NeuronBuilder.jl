using ModelingToolkit, OrdinaryDiffEq, Plots, .NeuronBuilder, Unitful
import Unitful: mV, mS, cm, mm, nA, mA, µA, ms, nF, μM

const area = Area(0.628e-3cm^2) # Prinz/Liu 0.0628 mm^2
const Cₘ = Cm() # specific capacitance cₘ is a biological constant (around) 10 nF/mm^2

Prinz_conv = Cₘ / area
syn_conv_factor = 1e-3 / area^2 #this gives synaptic conductances in μS/mm^2 

ics = ICS(Dict(:V => -60.0mV, :Ca => 0.05μM))

reversals = Reversals(Dict(:ENa => 50.0mV, :EH => -20.0mV, :EK => -80.0mV, :ELeak => -50.0mV))

params = Params(Dict(:area => area, :Cₘ => Cₘ, :τCa => 200.0ms, :Ca∞ => 0.05µM, :Iapp => 0.0nA/mm^2))

compart = Soma(ics, merge(reversals, params))
neur = Neuron(compart, AB2_ch, :ABneuron)
#Prinz_conv = Prinz_conversion(compart)


# AB2_ch = [Prinz.Na(100 * Prinz_conv), Prinz.CaS(6 * Prinz_conv), Prinz.CaT(2.5 * Prinz_conv), Prinz.H(0.01 * Prinz_conv),
#     Prinz.Ka(50 * Prinz_conv), Prinz.KCa(5 * Prinz_conv), Prinz.Kdr(100 * Prinz_conv), Prinz.Leak(0.0)]


# prob = ODEProblem(neur.ODESystem, [], (0.0, 2500.0), [])
# @time sol_individual = solve(prob, Rosenbrock23())
