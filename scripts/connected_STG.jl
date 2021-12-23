using NeuronBuilder, ModelingToolkit, OrdinaryDiffEq, Plots

# Using parameters from Prinz (2004) Similar network activity from disparate circuit parameters
# These are specifically Figure 3e (and Table 2 for current values)

# Membrane ion channels

AB1_channels = [Prinz.NaV(100.0), Prinz.CaT(2.5), Prinz.CaS(6.0), Prinz.Ka(50.0), Prinz.KCa(5.0), Prinz.Kdr(100.0), Prinz.H(0.01), Prinz.Leak(0.00)]
LP1_channels = [Prinz.NaV(100.0), Prinz.CaT(0.0), Prinz.CaS(4.0), Prinz.Ka(20.0), Prinz.KCa(0.0), Prinz.Kdr(25.0), Prinz.H(0.05), Prinz.Leak(0.03)]
PY1_channels = [Prinz.NaV(100.0), Prinz.CaT(2.5), Prinz.CaS(2.0), Prinz.Ka(50.0), Prinz.KCa(0.0), Prinz.Kdr(125.0), Prinz.H(0.05), Prinz.Leak(0.01)]

# Synapses from pre-synaptic AB/PD -> post-synaptic target
conv_factor = 1e-9 / (0.628e-3 * 1e-3) # HACK: divide original nS by Prinz model area and convert to mS

ABLP_chol = Chol(30.0 * conv_factor)
ABPY_chol = Chol(3.0 * conv_factor)
ABLP_glut = Glut(30.0 * conv_factor)
ABPY_glut = Glut(10.0 * conv_factor)

LPAB_glut = Glut(30.0 * conv_factor)
LPPY_glut = Glut(1.0 * conv_factor)

PYLP_glut = Glut(30.0 * conv_factor)

comp1 = deepcopy(comp)
comp1.hooks = 1
comp3 = deepcopy(comp)
comp3.hooks = 3
AB1 = build_neuron(comp1, AB1_channels; name = :AB1) #num_inputs = 1
PY1 = build_neuron(comp3, PY1_channels; name = :PY1) #num_inputs = 3
LP1 = build_neuron(comp3, LP1_channels; name = :LP1) #num_inputs = 3

grp = build_group([AB1, PY1, LP1]; name = :STG)

grp = add_connection(grp, AB1, LP1, ABLP_chol; i = 1)
grp = add_connection(grp, AB1, LP1, ABLP_glut; i = 2)
grp = add_connection(grp, PY1, LP1, PYLP_glut; i = 3)

grp = add_connection(grp, LP1, AB1, LPAB_glut; i = 1)

grp = add_connection(grp, LP1, PY1, LPPY_glut; i = 1)
grp = add_connection(grp, AB1, PY1, ABPY_glut; i = 2)
grp = add_connection(grp, AB1, PY1, ABPY_chol; i = 3)


tspan = (0.0, 5000.0)
final_sys = structural_simplify(grp)

prob = ODEProblem(final_sys, [], tspan, [])

sol = solve(prob, ImplicitEuler(), dt = 0.025)

plot(sol)


