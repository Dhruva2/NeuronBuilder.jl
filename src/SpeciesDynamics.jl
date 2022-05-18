"""

abstract type SpeciesDynamics{F<:Species} end

Each SpeciesDynamics struct needs a functor: 
    (s::SomeDynamics)(n::Neuron, variable, currents)

"""


get_parameters(::SpeciesDynamics) = nothing
get_states(::SpeciesDynamics) = nothing
default_params(::SpeciesDynamics, a, b, c) = Dict{Num,Float64}()
default_states(::SpeciesDynamics, a, b, c) = Dict{Num,Float64}()


struct BasicVoltageDynamics <: SpeciesDynamics{Voltage} end

function (b::BasicVoltageDynamics)(n::Neuron, vars, varnames, flux)
    Cₘ, = get_parameters(b)
    V = vars[findfirst(x -> x == Voltage, varnames)]
    return D(V) ~ (1 / Cₘ) * (flux) #non-standard convention, sum of fluxes already has negative sign because of the (E-V) in currents
end

get_parameters(::BasicVoltageDynamics) = @parameters Cₘ
default_params(v::BasicVoltageDynamics, n::Neuron, vars, varnames) = Dict(get_parameters(v)... => capacitance(n.geometry))