"""
Running the STG network 
- For units to check out, synapses also get a conversion factor.
- The names of the variables in the system can be listed with `@show stg.states`. 
- You can plot variables by their names with plot(sol;  vars = [STG.AB₊V]). 
"""

using ModelingToolkit, OrdinaryDiffEq, Plots

const Area = 0.0628 # Prinz/Liu 0.0628 mm2
const Cm = 10.0 # specific capacitance cₘ is a biological constant (around) 10 nF/mm^2

Prinz_conv = Cm / Area
synaptic_conv = 1e-3 / Area^2 #this gives μS/mm^2 

τCa = 200.0
Ca∞ = 0.05

defaults = Dict(Voltage => BasicVoltageDynamics(),
    Calcium => Prinz.CalciumDynamics(τCa, Ca∞, 14.96 * 0.0628),
    Reversal{Calcium} => Prinz.CaReversalDynamics())

somatic_parameters = Dict(
    Reversal{Sodium} => 50.0,
    Reversal{Potassium} => -80.0,
    Reversal{Leak} => -50.0,
    Reversal{Proton} => -20.0,
    Voltage => -60.0,
    Calcium => 0.05,
    Reversal{Calcium} => 0.0)

#Using parameters from Prinz (2004) Similar network activity from disparate circuit parameters
#these are specifically in Table 2: AB2, LP4, PY1 for figure 3e
function AB_Neuron(num_connections)
    AB2_ch = [
        Prinz.Na(100.0 * Prinz_conv),
        Prinz.CaS(6.0 * Prinz_conv),
        Prinz.CaT(2.5 * Prinz_conv),
        Prinz.H(0.01 * Prinz_conv),
        Prinz.Ka(50.0 * Prinz_conv),
        Prinz.KCa(5.0 * Prinz_conv),
        Prinz.Kdr(100.0 * Prinz_conv),
        Prinz.leak(0.0 * Prinz_conv)
    ]
    b = BasicNeuron(NoGeometry(Cm), defaults, somatic_parameters, AB2_ch, :AB)
    return b(incoming_connections=num_connections)
end

function PY_Neuron(num_connections)
    PY1_ch = [Prinz.Na(100.0 * Prinz_conv),
        Prinz.CaS(2.0 * Prinz_conv),
        Prinz.CaT(2.4 * Prinz_conv),
        Prinz.H(0.05 * Prinz_conv),
        Prinz.Ka(50.0 * Prinz_conv),
        Prinz.KCa(0.0 * Prinz_conv),
        Prinz.Kdr(125.0 * Prinz_conv),
        Prinz.leak(0.01 * Prinz_conv)
    ]
    somatic_parameters[Voltage] = -55.0
    b = BasicNeuron(NoGeometry(Cm), defaults, somatic_parameters, PY1_ch, :PY)
    return b(incoming_connections=num_connections)
end

function LP_Neuron(num_connections)
    LP4_ch = [Prinz.Na(100.0 * Prinz_conv),
        Prinz.CaS(4.0 * Prinz_conv),
        Prinz.CaT(0.0 * Prinz_conv),
        Prinz.H(0.05 * Prinz_conv),
        Prinz.Ka(20.0 * Prinz_conv),
        Prinz.KCa(0.0 * Prinz_conv),
        Prinz.Kdr(25.0 * Prinz_conv),
        Prinz.leak(0.03 * Prinz_conv)
    ]
    somatic_parameters[Voltage] = -65.0
    b = BasicNeuron(NoGeometry(Cm), defaults, somatic_parameters, LP4_ch, :LP)
    return b(incoming_connections=num_connections)
end

function nodes(i)
    (i == 1) && return n -> AB_Neuron(n)
    (i == 2) && return n -> PY_Neuron(n)
    (i == 3) && return n -> LP_Neuron(n)
end

intoAB = [(EmptyConnection(),)
    (EmptyConnection(),) #hi
    (Glut(30.0 * synaptic_conv),)
]

intoPY = [
    (Chol(3.0 * synaptic_conv), Glut(10.0 * synaptic_conv))
    (EmptyConnection(),)
    (Glut(1.0 * synaptic_conv),)
]

intoLP = [
    (Chol(30.0 * synaptic_conv), Glut(30.0 * synaptic_conv))
    (Glut(30.0 * synaptic_conv),)
    (EmptyConnection(),)
]

Adjacency_matrix = cat(intoAB, intoPY, intoLP; dims=2)

function edges(i, j)
    return Adjacency_matrix[i, j]
end

tspan = (0.0, 10000.0)
stg = @time build_network(nodes, edges, 1:3, name=:stg)
prob = ODEProblem(stg, [], tspan, []; jac=true)

@time sol = solve(prob, AutoTsit5(Rosenbrock23()))

plot(sol, xlims=(5000, 10000), fg_legend=:transparent; vars=[stg.AB₊V, stg.LP₊V, stg.PY₊V, stg.PY₊CaS₊ICa], layout=(4, 1))


"""
Another way of adding connections without having to specify a connectivity matrix (no empty connections)
"""
AB = AB_Neuron(1)
LP = LP_Neuron(3)
PY = PY_Neuron(3)
neurons = [AB, PY, LP]

ABLP_chol = directed_synapse(AB, LP, Chol(30.0 * synaptic_conv))
ABPY_chol = directed_synapse(AB, PY, Chol(3.0 * synaptic_conv))
ABLP_glut = directed_synapse(AB, LP, Glut(30.0 * synaptic_conv))
ABPY_glut = directed_synapse(AB, PY, Glut(10.0 * synaptic_conv))

LPAB_glut = directed_synapse(LP, AB, Glut(30.0 * synaptic_conv))
LPPY_glut = directed_synapse(LP, PY, Glut(1.0 * synaptic_conv))

PYLP_glut = directed_synapse(PY, LP, Glut(30.0 * synaptic_conv))

connections = [LPAB_glut, ABPY_chol, ABPY_glut, LPPY_glut, ABLP_chol, ABLP_glut, PYLP_glut]

stg_2 = @time add_all_connections(neurons, connections)
prob_2 = ODEProblem(stg_2, [], tspan, []; jac=true)

@time sol_2 = solve(prob_2, AutoTsit5(Rosenbrock23()))

plot(sol_2, xlims=(5000, 10000), fg_legend = :transparent; idxs=[stg_2.AB₊V, stg_2.LP₊V, stg_2.PY₊V, stg_2.PY₊CaS₊ICa], layout=(4, 1))


