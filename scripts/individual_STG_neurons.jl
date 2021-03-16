using NeuronBuilder, ModelingToolkit, OrdinaryDiffEq, Plots


AB1_channels = [NaV(100.), CaS(3.), CaT(1.3), H(0.5), Ka(5.), KCa(10.), Kdr(20.), Leak(0.01)]

PY1_channels = [NaV(100.), CaS(1.), CaT(1.3), H(2.5), Ka(5.), KCa(0.), Kdr(25.), Leak(0.01)]


LP1_channels = [NaV(100.), CaS(2.), CaT(0.), H(2.5), Ka(2.), KCa(0.), Kdr(5.), Leak(0.03)]


PY1 = build_neuron(PY1_channels; name = :PY1)
LP1 = build_neuron(PY1_channels; name = :LP1)
AB1 = build_neuron(AB1_channels; name=:AB1)

final_PY1 = structural_simplify(PY1)
probPY1 = ODEProblem(final_PY1, 
                collect(ModelingToolkit.get_default_u0(final_PY1)),
                (0.,500.),
                collect(ModelingToolkit.get_default_p(final_PY1)))



final_LP1 = structural_simplify(LP1)
probLP1 = ODEProblem(final_LP1, 
                collect(ModelingToolkit.get_default_u0(final_LP1)),
                (0.,500.),
                collect(ModelingToolkit.get_default_p(final_LP1)))

final_AB1 = structural_simplify(AB1)
probAB1 = ODEProblem(final_AB1, 
                collect(ModelingToolkit.get_default_u0(final_AB1)),
                (0.,500.),
                collect(ModelingToolkit.get_default_p(final_AB1)))