# first run a single neuron to load everything with `include("individual_neurons.jl")`

using NeuronBuilder, ModelingToolkit, OrdinaryDiffEq, Plots

#mRNAs are regulated according to the equations in OLeary et. al. 2014 (https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4109293/)
taus = Dict(:τNa => 1e6, :τCaS => 1e6, :τCaT => 1e6, :τKCa => 1e6, :τKa => 1e6, :τKdr => 1e6, :τH => 1e6)

#fig 3.e
AB2_ch = Regulated.([Prinz.Na(100.0 * Prinz_conv), Prinz.CaS(6.0 * Prinz_conv), Prinz.CaT(2.5 * Prinz_conv), Prinz.H(0.01 * Prinz_conv), Prinz.Ka(50.0 * Prinz_conv), Prinz.KCa(5.0 * Prinz_conv), Prinz.Kdr(100.0 * Prinz_conv), Prinz.Leak(0.0 * Prinz_conv)])
LP4_ch = Regulated.([Prinz.Na(100.0 * Prinz_conv), Prinz.CaS(4.0 * Prinz_conv), Prinz.CaT(0.0 * Prinz_conv), Prinz.H(0.05 * Prinz_conv), Prinz.Ka(20.0 * Prinz_conv), Prinz.KCa(0.0 * Prinz_conv), Prinz.Kdr(25.0 * Prinz_conv), Prinz.Leak(0.03 * Prinz_conv)])
PY1_ch = Regulated.([Prinz.Na(100.0 * Prinz_conv), Prinz.CaS(2.0 * Prinz_conv), Prinz.CaT(2.4 * Prinz_conv), Prinz.H(0.05 * Prinz_conv), Prinz.Ka(50.0 * Prinz_conv), Prinz.KCa(0.0 * Prinz_conv), Prinz.Kdr(125.0 * Prinz_conv), Prinz.Leak(0.01 * Prinz_conv)])

reversals = Dict(:ENa => 50.0, :EH => -20.0, :EK => -80.0, :ELeak => -50.0)
comp1 = Soma(Dict(:V => -60.0, :Ca => 0.05), merge(reversals, params, taus), 1)
reversals = Dict(:ENa => 50.0, :EH => -20.0, :EK => -80.0, :ELeak => -50.0)
comp2 = Soma(Dict(:V => -55.0, :Ca => 0.05), merge(reversals, params, taus), 3)
reversals = Dict(:ENa => 50.0, :EH => -20.0, :EK => -80.0, :ELeak => -50.0)
comp3 = Soma(Dict(:V => -65.0, :Ca => 0.05), merge(reversals, params, taus), 3)

AB = build_neuron(comp1, AB2_ch; name=:AB) 
PY = build_neuron(comp2, PY1_ch; name=:PY) 
LP = build_neuron(comp3, LP4_ch; name=:LP) 

grp = build_group([AB, PY, LP]; name=:STG)

conv_factor = syn_conv_factor(comp1)
#Synapses from fig 3.e
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
grp_final = structural_simplify(grp)
prob = ODEProblem(grp_final, [], tspan, [])

@time sol = solve(prob, AutoTsit5(Rosenbrock23()))

plot(sol, legend=false, xlims=(5000, 10000), ylims=(-80, 70))

# every regulated neuron has 27 variables, we get only volatge with 1, 28, 55
plot(sol.t, sol(sol.t, idxs=1), ylims=(-80, 70), label="AB")
plot!(sol.t, sol(sol.t, idxs=28), ylims=(-80, 70), label="PY")
plot!(sol.t, sol(sol.t, idxs=55), ylims=(-80, 70), label="LP")