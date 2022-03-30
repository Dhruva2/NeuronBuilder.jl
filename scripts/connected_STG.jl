```
Running the STG network 
- First run a single neuron to load everything with `include("individual_neurons.jl")` in the REPL
- Synapses also get a conversion factor which depends on the geometry of the somatic compartment.
- The group of neurons and their connections is an `ODESystem` on which we do `structural_simplify(grp)` to eliminate redundant variables. 
- The names of the variables in the system can be listed with `@show grp.states`. The problem that then takes this system is a regular `ODEProblem`.
- You can plot variables by their names with plot(sol;  vars = [STG.AB₊V]). 
```
using NeuronBuilder, ModelingToolkit, OrdinaryDiffEq, Plots

#Using parameters from Prinz (2004) Similar network activity from disparate circuit parameters
#these are specifically in Table 2: AB2, LP4, PY1 for figure 3e

#fig 3.e
AB2_ch = [Prinz.Na(100.0 * Prinz_conv), Prinz.CaS(6.0 * Prinz_conv), Prinz.CaT(2.5 * Prinz_conv), Prinz.H(0.01 * Prinz_conv), Prinz.Ka(50.0 * Prinz_conv), Prinz.KCa(5.0 * Prinz_conv), Prinz.Kdr(100.0 * Prinz_conv), Prinz.Leak(0.0 * Prinz_conv)]
LP4_ch = [Prinz.Na(100.0 * Prinz_conv), Prinz.CaS(4.0 * Prinz_conv), Prinz.CaT(0.0 * Prinz_conv), Prinz.H(0.05 * Prinz_conv), Prinz.Ka(20.0 * Prinz_conv), Prinz.KCa(0.0 * Prinz_conv), Prinz.Kdr(25.0 * Prinz_conv), Prinz.Leak(0.03 * Prinz_conv)]
PY1_ch = [Prinz.Na(100.0 * Prinz_conv), Prinz.CaS(2.0 * Prinz_conv), Prinz.CaT(2.4 * Prinz_conv), Prinz.H(0.05 * Prinz_conv), Prinz.Ka(50.0 * Prinz_conv), Prinz.KCa(0.0 * Prinz_conv), Prinz.Kdr(125.0 * Prinz_conv), Prinz.Leak(0.01 * Prinz_conv)]

reversals = Dict(:ENa => 50.0, :EH => -20.0, :EK => -80.0, :ELeak => -50.0)
comp1 = Soma(Dict(:V => -60.0, :Ca => 0.05), merge(reversals, params);num_connection= 1) 
reversals = Dict(:ENa => 50.0, :EH => -20.0, :EK => -80.0, :ELeak => -50.0)
comp2 = Soma(Dict(:V => -55.0, :Ca => 0.05), merge(reversals, params);num_connection=  3) 
reversals = Dict(:ENa => 50.0, :EH => -20.0, :EK => -80.0, :ELeak => -50.0)
comp3 = Soma(Dict(:V => -65.0, :Ca => 0.05), merge(reversals, params);num_connection=  3) 

AB = build_neuron(comp1, AB2_ch; name=:AB)
PY = build_neuron(comp2, PY1_ch; name=:PY)
LP = build_neuron(comp3, LP4_ch; name=:LP)

grp = build_group([AB, PY, LP]; name=:STG)

conv_factor = syn_conv_factor(comp1)
#Synapses 
ABLP_chol = Chol(30.0 * conv_factor)
ABPY_chol = Chol(3.0 * conv_factor)
ABLP_glut = Glut(30.0 * conv_factor)
ABPY_glut = Glut(10.0 * conv_factor)

LPAB_glut = Glut(30.0 * conv_factor)
LPPY_glut = Glut(1.0 * conv_factor)

PYLP_glut = Glut(30.0 * conv_factor)

grp = add_connection(grp, AB, LP, ABLP_chol; i=1)
grp = add_connection(grp, AB, LP, ABLP_glut; i=2)
grp = add_connection(grp, PY, LP, PYLP_glut; i=3)

grp = add_connection(grp, LP, AB, LPAB_glut; i=1)

grp = add_connection(grp, LP, PY, LPPY_glut; i=1)
grp = add_connection(grp, AB, PY, ABPY_glut; i=2)
grp = add_connection(grp, AB, PY, ABPY_chol; i=3)

tspan = (0.0, 10000.0)
STG = structural_simplify(grp)
prob = ODEProblem(STG, [], tspan, [])

@time sol = solve(prob, AutoTsit5(Rosenbrock23()))

plot(sol, xlims=(5000, 10000), ylims=(-80, 70); vars = [STG.AB₊V, STG.LP₊V, STG.PY₊V ])
