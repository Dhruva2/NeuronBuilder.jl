using NeuronBuilder, ModelingToolkit, OrdinaryDiffEq, Plots

# Using parameters from Prinz (2004) Similar network activity from disparate circuit parameters
# These are specifically Figure 3e (and Table 2 for current values)

# Membrane ion channels
AB1_channels = [NaV(100.), CaT(2.5), CaS(6.), Ka(50.), KCa(5.), Kdr(100.), H(0.01), Leak(0.00)]
LP1_channels = [NaV(100.), CaT(0.0), CaS(4.), Ka(20.), KCa(0.), Kdr(25.), H(0.05), Leak(0.03)]
PY1_channels = [NaV(100.), CaT(2.5), CaS(2.), Ka(50.), KCa(0.), Kdr(125.), H(0.05), Leak(0.01)]

# Synapses from pre-synaptic AB/PD -> post-synaptic target
conv_factor = 1e-9/(0.628e-3*1e-3) # HACK: divide original nS by Prinz model area and convert to mS

ABLP_chol = Chol(30.0*conv_factor)  
ABPY_chol = Chol(3.0*conv_factor) 
ABLP_glut = Glut(30.0*conv_factor)  
ABPY_glut = Glut(10.0*conv_factor)  

LPAB_glut = Glut(30.0*conv_factor)  
LPPY_glut = Glut(1.0*conv_factor) 

PYLP_glut = Glut(30.0*conv_factor)  

AB1 = build_neuron(AB1_channels; name = :AB1, num_inputs = 1)
PY1 = build_neuron(PY1_channels; name = :PY1, num_inputs = 3)
LP1 = build_neuron(LP1_channels; name = :LP1, num_inputs = 3)

grp = build_group([AB1, PY1, LP1]; name=:STG)

grp = add_connection(grp, AB1, LP1, ABLP_chol; i=1)
grp = add_connection(grp, AB1, LP1, ABLP_glut; i=2)
grp = add_connection(grp, PY1, LP1, PYLP_glut; i=3)

grp = add_connection(grp, LP1, AB1, LPAB_glut; i=1)

grp = add_connection(grp, LP1, PY1, LPPY_glut; i=1)
grp = add_connection(grp, AB1, PY1, ABPY_glut; i=2)
grp = add_connection(grp, AB1, PY1, ABPY_chol; i=3)


tspan = (0., 5000.)
final_sys = structural_simplify(grp)

prob = ODEProblem(final_sys, [], tspan,[])

sol = solve(prob, ImplicitEuler(), dt = 0.025)

plot(sol)
