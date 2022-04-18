using NeuronBuilder, ModelingToolkit, OrdinaryDiffEq, Plots

const Area = 0.0628 # Prinz/Liu 0.0628 mm2
const Cm = 10.0 # specific capacitance cₘ is a biological constant (around) 10 nF/mm^2

Liu_conv = Cm

channels = [Liu.Na(700.0 * Liu_conv), Liu.CaS(4.0 * Liu_conv), Liu.CaT(2.0 * Liu_conv), Liu.Ka(50.0 * Liu_conv), Liu.KCa(40.0 * Liu_conv),
    Liu.Kdr(70.0 * Liu_conv), Liu.H(0.03 * Liu_conv), Liu.leak(0.01 * Liu_conv)]
τmRNAs = [1e6 for el in channels] #could be different
τg = 1e6
Ca_tgt = 5.0
Oreg_channels = OLeary_reg.(channels, τmRNAs, τg, Ca_tgt)
A = rand(8,8)
#Freg_channels = [Franci_reg.(ch, τmRNAs, τg, Ca_tgt, A, channels) for ch in channels]

τCa = 20.0
Ca∞ = 0.05
fxarea = 14.96 * 0.0628

defaults = Dict(Voltage => BasicVoltageDynamics(),
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

b = BasicNeuron(NoGeometry(Cm), defaults, somatic_parameters, Oreg_channels, :regul_Liu)

neur = b(0)

prob = ODEProblem(neur.sys, [], (0.0, 5000.0), [])

sol = solve(prob, Tsit5())

sys = neur.sys
plot(sol, vars=[sys.V, sys.Ca])
plot(sol, vars=[sys.Na_regul₊mRNA, sys.CaS_regul₊mRNA, sys.CaT_regul₊mRNA, sys.Ka_regul₊mRNA, sys.KCa_regul₊mRNA, sys.Kdr_regul₊mRNA, sys.H_regul₊mRNA])