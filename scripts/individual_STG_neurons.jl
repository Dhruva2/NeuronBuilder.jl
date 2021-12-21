using NeuronBuilder, OrdinaryDiffEq, ModelingToolkit
using Plots
using NeuronBuilder.Prinz

ics = Dict(:V => -60.0, :Ca => 0.05)
reversals = Dict(:ENa => 50.0, :EH => -20.0, :EK => -80.0, :ELeak => -50.0)
params = Dict(:Cₘ => 10.0, :τCa => 200.0, :Ca∞ => 0.05, :Ca_tgt => 50.0, :τg => 1.0e6)
comp = Soma(ics, merge(reversals, params), 0)


# ## xolotl versions divided by 10
AB1_channels = [Prinz.NaV(1000.0), Prinz.CaS(60.0), Prinz.CaT(25.0), Prinz.H(0.1), Prinz.Ka(500.0), Prinz.KCa(50.0), Prinz.Kdr(1000.0), Prinz.Leak(0.0)]

#LP_channels = [Prinz.NaV(1000.), Prinz.CaS(40.), Prinz.CaT(0.), Prinz.H(0.5), Prinz.Ka(200.), Prinz.KCa(0.), Prinz.Kdr(250.), Prinz.Leak(0.3)]

#PY1_channels = [Prinz.NaV(100.0), Prinz.CaS(2.0), Prinz.CaT(2.4), Prinz.H(0.05), Prinz.Ka(50.0),Prinz.KCa(0.0), Prinz.Kdr(125.0), Prinz.Leak(0.01)]

#PY1 = build_neuron(comp, PY1_channels, false; name = :PY1)
#LP = build_neuron(comp, LP1_channels, false; name = :LP1)
AB = build_neuron(comp, AB1_channels, false; name = :AB1)

final_AB = structural_simplify(AB)
probAB = ODEProblem(final_AB, [], (0.0, 3000.0), [])

#sol = solve(probAB1, ImplicitEuler(), dt = 0.025)
sol = solve(probAB, AutoTsit5(TRBDF2()))
plot(sol)
