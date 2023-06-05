using ModelingToolkit, OrdinaryDiffEq, Plots

const Cm = 10.0 # specific capacitance cₘ is a biological constant (around) 10 nF/mm^2
Liu_conv = Cm

channels = [Liu.Na(700.0 * Liu_conv), Liu.CaS(4.0 * Liu_conv), Liu.CaT(2.0 * Liu_conv), Liu.Ka(50.0 * Liu_conv), Liu.KCa(40.0 * Liu_conv),
    Liu.Kdr(70.0 * Liu_conv), Liu.H(0.03 * Liu_conv), Liu.leak(0.01 * Liu_conv)]

τmRNAs = repeat([5e4], 7)
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

b = BasicNeuron(NoGeometry(Cm), defaults, somatic_parameters, Oreg_channels, :regul_Liu)

neur = b()

prob = ODEProblem(neur, [], (0.0, 15000.0), [])

sol = @time solve(prob, AutoTsit5(Rosenbrock23()))


p1 = plot(sol, idxs=[neur.V, neur.Ca])
p2 = plot(sol, idxs=[neur.Na_regul₊mRNA, neur.CaS_regul₊mRNA, neur.CaT_regul₊mRNA, neur.Ka_regul₊mRNA, neur.KCa_regul₊mRNA, neur.Kdr_regul₊mRNA, neur.H_regul₊mRNA], title="mRNAs", titlefontsize=9)
p3 = plot(sol, idxs=[neur.Na_regul₊gNa, neur.CaS_regul₊gCa, neur.CaT_regul₊gCa, neur.Ka_regul₊gK, neur.KCa_regul₊gK, neur.Kdr_regul₊gK, neur.H_regul₊gH], title="conductances", titlefontsize=9)
plot(p1,p2,p3,layout=(3,1), legend=false)