using  ModelingToolkit, OrdinaryDiffEq, Plots

## xolotl versions
# AB1_channels = [NaV(1000.), CaS(60.), CaT(25.), H(0.1), Ka(500.), KCa(50.), Kdr(1000.), Leak(0.)]

# LP1_channels = [NaV(1000.), CaS(40.), CaT(0.), H(0.5), Ka(200.), KCa(0.), Kdr(250.), Leak(0.3)]

# PY1_channels = [NaV(1000.), CaS(20.), CaT(24.), H(0.5), Ka(500.), KCa(0.), Kdr(1250.), Leak(0.1)]


## xolotl versions divided by 10
AB1_channels = [NaV(100.), CaS(6.), CaT(2.5), H(0.01), Ka(50.), KCa(5.), Kdr(100.), Leak(0.)]

LP1_channels = [NaV(100.), CaS(4.), CaT(0.), H(0.05), Ka(20.), KCa(0.), Kdr(25.), Leak(0.03)]

PY1_channels = [NaV(100.), CaS(2.), CaT(2.4), H(0.05), Ka(50.), KCa(0.), Kdr(125.), Leak(0.01)]


## ab1 channels that actually works
AB1_channels = [NaV(100.), CaS(3.), CaT(1.3), H(0.5), Ka(5.), KCa(10.), Kdr(20.), Leak(0.01)]



PY1 = build_neuron(PY1_channels; name = :PY1)
LP1 = build_neuron(PY1_channels; name = :LP1)
AB1 = build_neuron(AB1_channels; name = :AB1)

final_PY1 = structural_simplify(PY1)
probPY1 = ODEProblem(final_PY1, [], (0., 3000.), [])



final_LP1 = structural_simplify(LP1)
probLP1 = ODEProblem(final_LP1,[],
                (0.,3000.), [])

final_AB1 = structural_simplify(AB1)

probAB1 = ODEProblem(final_AB1, [], (0., 3000.), [])

