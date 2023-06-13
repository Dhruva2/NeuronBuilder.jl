


using ModelingToolkit, OrdinaryDiffEq, Plots

Liu_conv = 10.0

liu_channels = [Liu.Na(700.0 * Liu_conv), Liu.CaS(4.0 * Liu_conv), Liu.CaT(2.0 * Liu_conv), Liu.Ka(50.0 * Liu_conv), Liu.KCa(40.0 * Liu_conv),
    Liu.Kdr(70.0 * Liu_conv), Liu.H(0.03 * Liu_conv), Liu.Leak(0.01 * Liu_conv)]


τCa = 20.0
Ca∞ = 0.05
fxarea = 14.96 * 0.0628


# synapses = [Liu.chol(4.0)]
synapses = [EmptySynapse()]

liu() = BasicNeuron(
    :liu,
    NoGeometry(Dict(Capacitance{Voltage}() => 10.0)),
    Dict(
        Voltage() => BasicComponents.basic_influx_dynamics(Voltage()),
        Calcium() => Liu.calcium_dynamics(; τCa=20.0, Ca∞=0.05, calc_multiplier=14.96 * 0.0628),
        Reversal{Calcium}() => Liu.calcium_reversal_dynamics()
    ),
    Dict(
        Reversal{Sodium}() => 50.0,
        Reversal{Potassium}() => -80.0,
        Reversal{Voltage}() => -50.0,
        Reversal{Proton}() => -20.0,
        Voltage() => -60.0,
        Calcium() => 0.05,
        Reversal{Calcium}() => 0.0),
    liu_channels,
    synapses
)

b = liu()
sys = b() |> structural_simplify
prob = ODEProblem(sys, [], (0.0, 5000.0), [])
sol = solve(prob, Tsit5())
plot(sol)
