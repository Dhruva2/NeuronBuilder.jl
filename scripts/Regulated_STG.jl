using NeuronBuilder, ModelingToolkit, OrdinaryDiffEq, Plots


const Area = 0.0628 # Prinz/Liu 0.0628 mm2
const Cm = 10.0 # specific capacitance cₘ is a biological constant (around) 10 nF/mm^2

Prinz_conv = Cm / Area
synaptic_conv = 1e-3 / Area^2 #this gives μS/mm^2 

params = Dict(:area => Area, :Cₘ => Cm, :τCa => 200.0, :Ca∞ => 0.05, :Ca_tgt => 200.0, :τg => 1e6, :Iapp => 0.0)

reversals = Dict(:ENa => 50.0, :EH => -20.0, :EK => -80.0, :ELeak => -50.0)


#mRNAs are regulated according to the equations in OLeary et. al. 2014 (https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4109293/)
taus = Dict(:τNa => 1e6, :τCaS => 1e6, :τCaT => 1e6, :τKCa => 1e6, :τKa => 1e6, :τKdr => 1e6, :τH => 1e6)

#fig 3.e
AB2_ch = Regulated.([Prinz.Na(100.0 * Prinz_conv), Prinz.CaS(6.0 * Prinz_conv), Prinz.CaT(2.5 * Prinz_conv), Prinz.H(0.01 * Prinz_conv), Prinz.Ka(50.0 * Prinz_conv), Prinz.KCa(5.0 * Prinz_conv), Prinz.Kdr(100.0 * Prinz_conv), Prinz.Leak(0.0 * Prinz_conv)])
LP4_ch = Regulated.([Prinz.Na(100.0 * Prinz_conv), Prinz.CaS(4.0 * Prinz_conv), Prinz.CaT(0.0 * Prinz_conv), Prinz.H(0.05 * Prinz_conv), Prinz.Ka(20.0 * Prinz_conv), Prinz.KCa(0.0 * Prinz_conv), Prinz.Kdr(25.0 * Prinz_conv), Prinz.Leak(0.03 * Prinz_conv)])
PY1_ch = Regulated.([Prinz.Na(100.0 * Prinz_conv), Prinz.CaS(2.0 * Prinz_conv), Prinz.CaT(2.4 * Prinz_conv), Prinz.H(0.05 * Prinz_conv), Prinz.Ka(50.0 * Prinz_conv), Prinz.KCa(0.0 * Prinz_conv), Prinz.Kdr(125.0 * Prinz_conv), Prinz.Leak(0.01 * Prinz_conv)])


comp1 = Soma(Dict(:V => -60.0, :Ca => 0.05), merge(reversals, params, taus))
comp2 = Soma(Dict(:V => -55.0, :Ca => 0.05), merge(reversals, params, taus))
comp3 = Soma(Dict(:V => -65.0, :Ca => 0.05), merge(reversals, params, taus))

AB = Neuron(comp1, AB2_ch, 1, :AB)
PY = Neuron(comp2, PY1_ch, 3, :PY)
LP = Neuron(comp3, LP4_ch, 3, :LP)

grp = build_group([AB, PY, LP]; name=:STG)


#Synapses from fig 3.e
ABLP_chol = Chol(30.0 * synaptic_conv)
ABPY_chol = Chol(3.0 * synaptic_conv)
ABLP_glut = Glut(30.0 * synaptic_conv)
ABPY_glut = Glut(10.0 * synaptic_conv)

LPAB_glut = Glut(30.0 * synaptic_conv)
LPPY_glut = Glut(1.0 * synaptic_conv)

PYLP_glut = Glut(30.0 * synaptic_conv)

grp = add_connection(grp, AB, LP, ABLP_chol; i=1)
grp = add_connection(grp, AB, LP, ABLP_glut; i=2)
grp = add_connection(grp, PY, LP, PYLP_glut; i=3)

grp = add_connection(grp, LP, AB, LPAB_glut; i=1)

grp = add_connection(grp, LP, PY, LPPY_glut; i=1)
grp = add_connection(grp, AB, PY, ABPY_glut; i=2)
grp = add_connection(grp, AB, PY, ABPY_chol; i=3)

tspan = (0.0, 10000.0)
STG_R = structural_simplify(grp)
prob = ODEProblem(STG_R, [], tspan, [])

@time sol = solve(prob, AutoTsit5(Rosenbrock23()))

plot(sol, xlims=(0, 10000), ylims=(-80, 70); vars=[STG_R.AB₊V, STG_R.LP₊V, STG_R.PY₊V], layout=(3, 1))