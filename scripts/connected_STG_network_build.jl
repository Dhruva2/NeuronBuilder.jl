using NeuronBuilder, ModelingToolkit, OrdinaryDiffEq, Plots



const Area = 0.0628 # Prinz/Liu 0.0628 mm2
const Cm = 10.0 # specific capacitance cₘ is a biological constant (around) 10 nF/mm^2

Prinz_conv = Cm / Area
synaptic_conv = 1e-3 / Area^2 #this gives μS/mm^2 

params = Dict(:area => Area, :Cₘ => Cm, :τCa => 200.0, :Ca∞ => 0.05, :Ca_tgt => 200.0, :τg => 1e6, :Iapp => 0.0)

reversals = Dict(:ENa => 50.0, :EH => -20.0, :EK => -80.0, :ELeak => -50.0)


function AB_Neuron(num_connections)
    channels = [Prinz.Na(100.0 * Prinz_conv), Prinz.CaS(6.0 * Prinz_conv), Prinz.CaT(2.5 * Prinz_conv), Prinz.H(0.01 * Prinz_conv), Prinz.Ka(50.0 * Prinz_conv), Prinz.KCa(5.0 * Prinz_conv), Prinz.Kdr(100.0 * Prinz_conv), Prinz.Leak(0.0 * Prinz_conv)]
    comp = Soma(Dict(:V => -60.0, :Ca => 0.05), merge(reversals, params))
    return Neuron(comp, channels, num_connections, :AB)
end

function PY_Neuron(num_connections)
    channels = [Prinz.Na(100.0 * Prinz_conv), Prinz.CaS(2.0 * Prinz_conv), Prinz.CaT(2.4 * Prinz_conv), Prinz.H(0.05 * Prinz_conv), Prinz.Ka(50.0 * Prinz_conv), Prinz.KCa(0.0 * Prinz_conv), Prinz.Kdr(125.0 * Prinz_conv), Prinz.Leak(0.01 * Prinz_conv)]

    comp = Soma(Dict(:V => -55.0, :Ca => 0.05), merge(reversals, params))
    return Neuron(comp, channels, 3, :PY)
end

function LP_Neuron(num_connections)
    channels = [Prinz.Na(100.0 * Prinz_conv), Prinz.CaS(4.0 * Prinz_conv), Prinz.CaT(0.0 * Prinz_conv), Prinz.H(0.05 * Prinz_conv), Prinz.Ka(20.0 * Prinz_conv), Prinz.KCa(0.0 * Prinz_conv), Prinz.Kdr(25.0 * Prinz_conv), Prinz.Leak(0.03 * Prinz_conv)]
    comp = Soma(Dict(:V => -65.0, :Ca => 0.05), merge(reversals, params))
    return Neuron(comp, channels, 3, :LP)
end

# or return x -> Neuron(num_connections)
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
stg = build_network(nodes, edges, 3)
prob = ODEProblem(stg, [], tspan, [])

@time sol = solve(prob, AutoTsit5(Rosenbrock23()))

plot(sol, xlims=(5000, 10000), ylims=(-80, 70); vars=[stg.AB₊V, stg.LP₊V, stg.PY₊V], layout=(3, 1))
