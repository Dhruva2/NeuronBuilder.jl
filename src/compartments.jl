"""
To deal with calcium issue:

Take tau_ca out of the calcium currents themselves, and put them in the soma, like a capacitance


What do compartments need?
Reversal potentials
Capacitance
Hooks
Time constants (Ca + V)

"""

abstract type Compartment end

struct Soma{S} <: Compartment
    capacitance::S 
    reversals::Dict{Symbol, Float64}
    τCa::S
    hooks::Int64
end


