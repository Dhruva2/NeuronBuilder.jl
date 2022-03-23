using NeuronBuilder, OrdinaryDiffEq, ModelingToolkit, Plots

const Area = 0.0628 # Prinz/Liu 0.0628 mm2
const Cm = 10 # specific capacitance cₘ is a biological constant (around) 10 nF/mm^2

ics = Dict(:V => -60.0, :Ca => 0.05)
reversals = Dict(:ENa => 50.0, :EH => -20.0, :EK => -80.0, :ELeak => -50.0)

# [Ca]∞ units are uM. 
params = Dict(:area => Area, :Cₘ => Cm, :τCa => 200.0, :Ca∞ => 0.05, :Ca_tgt => 2.0, :τg => 1.0e6, :Iapp => 0.0)
taus = Dict(:τNa => 1e6, :τCaS => 1e6, :τCaT => 1e6, :τKCa => 1e6, :τKa => 1e6, :τKdr => 1e6, :τH => 1e6)
compart = Soma(ics, merge(reversals, params, taus))

Prinz_conv = Prinz_conversion(compart)

# pick a set of conductances from the file gvalues_papers.jl 
AB2_ch = Regulated.([Prinz.Na(100 * Prinz_conv), Prinz.CaS(6 * Prinz_conv), Prinz.CaT(2.5 * Prinz_conv), Prinz.H(0.01 * Prinz_conv),
    Prinz.Ka(50 * Prinz_conv), Prinz.KCa(5 * Prinz_conv), Prinz.Kdr(100 * Prinz_conv), Prinz.Leak(0.0)])

neur = build_neuron(compart, AB2_ch; name=:neuron)

prob = ODEProblem(neur, [], (0.0, 3000.0), [])
@time sol_individual = solve(prob, Rosenbrock23())

function my_plot(sol, name)
    p = plot(legend=:outertopright, title="$name")
    plot!(sol.t, sol(sol.t, idxs=1), linewidth=1, label="V", ylabel="mV")
    plot!(sol.t, sol(sol.t, idxs=2), linewidth=0.7, label="Ca [μM]")
    plot!(xlabel="Time (ms)")
    #savefig("$name.pdf")
end

my_plot(sol_individual, "AB Prinz")