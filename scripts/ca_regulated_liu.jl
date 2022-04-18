using NeuronBuilder, ModelingToolkit, OrdinaryDiffEq, Plots

const Area = 0.0628 # Prinz/Liu 0.0628 mm2
const Cm = 10.0 # specific capacitance cₘ is a biological constant (around) 10 nF/mm^2

Liu_conv = Cm

channels = [Liu.Na(700.0 * Liu_conv), Liu.CaS(4.0 * Liu_conv), Liu.CaT(2.0 * Liu_conv), Liu.Ka(50.0 * Liu_conv), Liu.KCa(40.0 * Liu_conv),
    Liu.Kdr(70.0 * Liu_conv), Liu.H(0.03 * Liu_conv), Liu.leak(0.01 * Liu_conv)]

τmRNAs = [5e4 5e4 5e4 5e4 5e4 5e4 5e4]
τg = 1e4
Ca_tgt = 5.0 

#average calcium at the start is 10, Ca_tgt is half as much so conductances have to decrease
reg_chs = [OLeary_reg(channel,taumrna,τg, Ca_tgt) for (channel, taumrna) in zip(channels[1:end-1],τmRNAs)]
Oreg_channels = [reg_chs... , channels[end]]

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

b = BasicNeuron(NoGeometry(Cm), defaults, somatic_parameters, Oreg_channels, :Oregul_Liu)

neur = b(0)

prob = ODEProblem(neur.sys, [], (0.0, 15000.0), [])

sol = @time solve(prob, AutoTsit5(Rosenbrock23()))


sys = neur.sys

p1 = plot(sol, vars=[sys.V, sys.Ca])
p2 = plot(sol, vars=[sys.Na_regul₊mRNA, sys.CaS_regul₊mRNA, sys.CaT_regul₊mRNA, sys.Ka_regul₊mRNA, sys.KCa_regul₊mRNA, sys.Kdr_regul₊mRNA, sys.H_regul₊mRNA], title="mRNAs", titlefontsize=9)
p3 = plot(sol, vars=[sys.Na_regul₊gNa, sys.CaS_regul₊gCa, sys.CaT_regul₊gCa, sys.Ka_regul₊gK, sys.KCa_regul₊gK, sys.Kdr_regul₊gK, sys.H_regul₊gH], title="conductances", titlefontsize=9)
plot(p1,p2,p3,layout=(3,1), legend=false)