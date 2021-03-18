using ModelingToolkit, OrdinaryDiffEq, Plots


# AB1_channels = [NaV(100.), CaS(3.), CaT(1.3), H(0.5), Ka(5.), KCa(10.), Kdr(20.), Leak(0.01)]
# PY1_channels = [NaV(100.), CaS(1.), CaT(1.3), H(2.5), Ka(5.), KCa(0.), Kdr(25.), Leak(0.01)]
# LP1_channels = [NaV(100.), CaS(2.), CaT(0.), H(2.5), Ka(2.), KCa(0.), Kdr(5.), Leak(0.03)]

AB1_channels = [NaV(100.), CaS(6.), CaT(2.5), H(0.01), Ka(50.), KCa(5.), Kdr(100.), Leak(0.)]
LP1_channels = [NaV(100.), CaS(4.), CaT(0.), H(0.05), Ka(20.), KCa(0.), Kdr(25.), Leak(0.03)]
PY1_channels = [NaV(100.), CaS(2.), CaT(2.4), H(0.05), Ka(50.), KCa(0.), Kdr(125.), Leak(0.01)]


## xolotl values
ABLP_chol = Chol(30.)
ABPY_chol = Chol(3.)
ABLP_glut = Glut(30.)
ABPY_glut = Glut(10.)
LPPY_glut = Glut(1.)
PYLP_glut = Glut(30.)
LPAB_glut = Glut(30.)

## xolotl values divided by 10
# ABLP_chol = Chol(3.)
# ABPY_chol = Chol(0.3)
# ABLP_glut = Glut(3.)
# ABPY_glut = Glut(1.)
# LPPY_glut = Glut(0.1)
# PYLP_glut = Glut(3.)
# LPAB_glut = Glut(3.)

ABLP_chol = Chol(0.03)
ABPY_chol = Chol(0.003)
ABLP_glut = Glut(0.03)
ABPY_glut = Glut(0.01)
LPPY_glut = Glut(0.001)
PYLP_glut = Glut(0.03)
LPAB_glut = Glut(0.03)


# ABLP_chol = Chol(3.)
# ABPY_chol = Chol(0.01)
# ABLP_glut = Glut(0.01)
# ABPY_glut = Glut(1.)
# LPPY_glut = Glut(0.1)
# PYLP_glut = Glut(0.1)
# LPAB_glut = Glut(10.)

AB1 = build_neuron(AB1_channels; name=:AB1, num_inputs = 1)
PY1 = build_neuron(PY1_channels; name = :PY1, num_inputs = 3)
LP1 = build_neuron(LP1_channels; name = :LP1, num_inputs = 3)

gr = build_group([AB1, PY1, LP1]; name=:STG)

gr = add_connection(gr, AB1, LP1, ABLP_chol; i=1)
gr = add_connection(gr, AB1, LP1, ABLP_glut; i=2)
gr = add_connection(gr, PY1, LP1, PYLP_glut; i=3)

gr = add_connection(gr, LP1, AB1, LPAB_glut; i=1)

gr = add_connection(gr, LP1, PY1, LPPY_glut; i=1)
gr = add_connection(gr, AB1, PY1, ABPY_glut; i=2)
gr = add_connection(gr, AB1, PY1, ABPY_chol; i=3)


tspan = (0., 1000.)
final_sys = structural_simplify(gr)

prob = ODEProblem(final_sys, 
                collect(ModelingToolkit.get_default_u0(final_sys)),
                tspan,
                collect(ModelingToolkit.get_default_p(final_sys)))
