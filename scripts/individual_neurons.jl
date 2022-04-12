using ModelingToolkit, OrdinaryDiffEq, Plots, NeuronBuilder

const Area = 0.0628 # Prinz/Liu 0.0628 mm2
const Cm = 10.0 # specific capacitance cₘ is a biological constant (around) 10 nF/mm^2

Liu_conv = Cm
Prinz_conv = Cm / Area
syn_conv_factor = 1e-3 / Area^2 #this gives μS/mm^2 

ics = Dict(:V => -60.0, :Ca => 0.05)
reversals = Dict(:ENa => 50.0, :EH => -20.0, :EK => -80.0, :ELeak => -50.0)

# [Ca]∞ units are uM. 
params = Dict(:area => Area, :Cₘ => Cm, :τCa => 200.0, :Ca∞ => 0.05, :Ca_tgt => 200.0, :τg => 1e6, :Iapp => 0.0)
compart = Soma(ics, merge(reversals, params))

DIC = [Liu.Na(700.0 * Liu_conv), Liu.CaS(4.0 * Liu_conv), Liu.CaT(2.0 * Liu_conv), Liu.Ka(50.0 * Liu_conv), Liu.KCa(40.0 * Liu_conv),
    Liu.Kdr(70.0 * Liu_conv), Liu.H(0.03 * Liu_conv), Liu.Leak(0.01 * Liu_conv)]
neur = Neuron(compart, DIC, 0, :Liuneuron)

prob = ODEProblem(neur.ODESystem, [], (0.0, 2500.0), [])
@time sol_individual = solve(prob, Rosenbrock23())

function my_plot(sol, name)
    p = plot(legend=:outertopright, title="$name")
    plot!(sol.t, sol(sol.t, idxs=1), linewidth=1, label="V", ylabel="mV")
    plot!(sol.t, sol(sol.t, idxs=2), linewidth=0.7, label="Ca [μM]")
    plot!(xlabel="Time (ms)")
    #savefig("$name.pdf")
end

my_plot(sol_individual, "Liu")