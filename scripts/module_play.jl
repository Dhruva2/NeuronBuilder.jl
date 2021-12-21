using NeuronBuilder, OrdinaryDiffEq, ModelingToolkit
using Plots

using NeuronBuilder.Liu

#main()
ics = Dict(:V => -60.0, :Ca => 0.05)
reversals = Dict(:ENa => 50.0, :EH => -20.0, :EK => -80.0, :ELeak => -50.0)
params = Dict(:Cₘ => 10.0, :τCa => 200.0, :Ca∞ => 0.05, :Ca_tgt => 50.0, :τg => 1e6)
#taus = Dict(:τNa => 1e6, :τCaS => 1e6, :τCaT => 1e6, :τKCa => 1e6, :τKa => 1e6, :τKdr => 1e6, :τH => 1e6)
comp = Soma(ics, merge(reversals, params), 0)

AB1_channels = [Liu.NaV(1831.0), Liu.CaS(27.0), Liu.CaT(23.0), Liu.H(10.1), Liu.Ka(246.0), Liu.KCa(980.0),
    Liu.Kdr(610.0), Liu.Leak(0.99)]

regulated = false
AB1 = build_neuron(comp, AB1_channels, regulated; name = :AB1)

final_AB1 = structural_simplify(AB1)
probAB1 = ODEProblem(final_AB1, [], (0.0, 3000.0), [])

#sol = solve(probAB1, ImplicitEuler(), dt = 0.025)
sol = solve(probAB1, AutoTsit5(TRBDF2()))
plot(sol)