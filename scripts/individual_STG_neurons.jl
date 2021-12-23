using NeuronBuilder, OrdinaryDiffEq, ModelingToolkit, Plots
using NeuronBuilder.Liu

ics = Dict(:V => -60.0, :Ca => 0.05)
reversals = Dict(:ENa => 50.0, :EH => -20.0, :EK => -80.0, :ELeak => -50.0)
params = Dict(:Cₘ => 10.0, :τCa => 200.0, :Ca∞ => 0.05, :Ca_tgt => 50.0, :τg => 1.0e6)
#taus = Dict(:τNa => 1e6, :τCaS => 1e6, :τCaT => 1e6, :τKCa => 1e6, :τKa => 1e6, :τKdr => 1e6, :τH => 1e6)
comp = Soma(ics, merge(reversals, params), 0)

##From MyModelMenagerie
#stg_liu = [Liu.NaV(100.0), Liu.CaS(3.0), Liu.CaT(1.3), Liu.Ka(5.0), Liu.KCa(10.0), Liu.Kdr(20.0), Liu.H(0.5), Liu.Leak(0.01)]

## parameters from Prinz 2003 https://journals.physiology.org/doi/pdf/10.1152/jn.00641.2003
#burster_channels = [Prinz.NaV(400.0), Prinz.CaS(8.0), Prinz.CaT(0.0), Prinz.H(0.04), Prinz.Ka(50.0), Prinz.KCa(20.0), Prinz.Kdr(50.0), Prinz.Leak(0.0)]

## parameters from O'Leary
AB1_channels = [Liu.NaV(1831.0), Liu.CaS(27.0), Liu.CaT(23.0), Liu.H(10.1), Liu.Ka(246.0), Liu.KCa(980.0), Liu.Kdr(610.0), Liu.Leak(0.99)]

AB = build_neuron(comp, AB1_channels; name = :AB1)
#burster = build_neuron(comp, stg_liu; name = :Burster)

final = structural_simplify(AB)
prob = ODEProblem(final, [], (0.0, 3000.0), [])

#sol = solve(probAB1, ImplicitEuler(), dt = 0.025)
sol = solve(prob, Tsit5())
plot(sol)
