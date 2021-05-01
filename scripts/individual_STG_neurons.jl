using NeuronBuilder, OrdinaryDiffEq, ModelingToolkit
using Plots; plotlyjs()

ics = Dict(:V => -60., :Ca => 0.05)
reversals = Dict(:ENa => 50., :EH =>-20., :EK =>-80., :ELeak =>-50.)
params = Dict(:Cₘ => 10., :τCa => 200., :Ca∞ => 0.05)
comp = Soma(ics, merge(reversals, params), 0)



# ## xolotl versions divided by 10
# AB1_channels = [NaV(100.), CaS(6.), CaT(2.5), H(0.01), Ka(50.), KCa(5.), Kdr(100.), Leak(0.)]

# LP1_channels = [NaV(100.), CaS(4.), CaT(0.), H(0.05), Ka(20.), KCa(0.), Kdr(25.), Leak(0.03)]

# PY1_channels = [NaV(100.), CaS(2.), CaT(2.4), H(0.05), Ka(50.), KCa(0.), Kdr(125.), Leak(0.01)]


## parameters from O'Leary
AB1_channels = [Liu_NaV(1831.), Liu_CaS(27.), Liu_CaT(23.), Liu_H(10.1), Liu_Ka(246.), Liu_KCa(980.), Liu_Kdr(610.), Liu_Leak(0.99)]


#PY1 = build_neuron(PY1_channels; name = :PY1)
#LP1 = build_neuron(PY1_channels; name = :LP1)
AB1 = build_neuron(comp, AB1_channels; name = :AB1)

final_AB1 = structural_simplify(AB1)
probAB1 = ODEProblem(final_AB1, [], (0., 3000.), [])

#sol = solve(probAB1, ImplicitEuler(), dt = 0.025)
sol = solve(probAB1, AutoTsit5(TRBDF2()))
plot(sol)
