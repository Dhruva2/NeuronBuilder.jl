using NeuronBuilder, OrdinaryDiffEq, ModelingToolkit, Plots

ics = Dict(:V => -60.0, :Ca => 0.05)
reversals = Dict(:ENa => 50.0, :EH => -20.0, :EK => -80.0, :ELeak => -5.0)
params = Dict(:Cₘ => 1.0, :τCa => 200.0, :Ca∞ => 0.05, :Ca_tgt => 50.0, :τg => 1.0e6, :Iapp => -0.5)
taus = Dict(:τNa => 1e6, :τCaS => 1e6, :τCaT => 1e6, :τKCa => 1e6, :τKa => 1e6, :τKdr => 1e6, :τH => 1e6)
compart = Soma(ics, merge(reversals, params, taus), 0)

##From MyModelMenagerie
#stg_liu = [Liu.NaV(100.0), Liu.CaS(3.0), Liu.CaT(1.3), Liu.Ka(5.0), Liu.KCa(10.0), Liu.Kdr(20.0), Liu.H(0.5), Liu.Leak(0.01)]
stg_liu = [Liu.NaV(55.0), Liu.CaS(2.0), Liu.CaT(1.4), Liu.Ka(52.0), Liu.KCa(30.0), Liu.Kdr(45.0), Liu.H(1.0), Liu.Leak(0.01)]
HH = [Liu.NaV(159.23), Liu.Leak(0.16), Liu.Kdr(47.77)]
tesis = [Liu.NaV(700.0), Liu.CaS(2.25), Liu.CaT(6.25), Liu.Ka(85.0), Liu.KCa(50.0), Liu.Kdr(90.0), Liu.H(0.03), Liu.Leak(0.01)]


## parameters from O'Leary
AB1_channels = [Liu.NaV(183.10), Liu.CaS(2.70), Liu.CaT(2.30), Liu.H(1.01), Liu.Ka(24.60), Liu.KCa(98.00), Liu.Kdr(61.00), Liu.Leak(0.099)]

AB = build_neuron(compart, HH; reg = true, name = :AB1)
#burster = build_neuron(comp, stg_liu; name = :Burster)

final = structural_simplify(AB)
prob = ODEProblem(final, [], (0.0, 3000.0), [])

#sol = solve(probAB1, ImplicitEuler(), dt = 0.025)
sol = solve(prob, Tsit5())
plot(sol)
