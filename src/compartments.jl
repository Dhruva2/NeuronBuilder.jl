"""
To deal with calcium issue:

Take tau_ca out of the calcium currents themselves, and put them in the soma, like a capacitance


What do compartments need?
Reversal potentials
Capacitance
Hooks
Time constants (Ca + V)
Dict(:V => -60., :Ca => 0.1),
Dict(:Cm => -60., :τCa => 10., :Cainf => 0.05, :ENa => -50.)
Dict(:Cₘ => -60., :τCa => 200., :Ca∞ => 0.05, :ENa => -50.)

[eNa, eh, eK, eleak
[50,-20,-80,-50
reversals = Dict(:ENa => 50., :EH =>-20., :EK =>80., :ELeak =>-50.)
ics = Dict(:V => -60., :Ca => 0.05)
params = Dict(:Cₘ => 10., :τCa => 200., :Ca∞ => 0.05)
"""

abstract type Compartment end

struct Soma <: Compartment 
    initial_states::Dict{Symbol, Float64}
    parameters::Dict{Symbol, Float64}
    hooks::Int64
end


