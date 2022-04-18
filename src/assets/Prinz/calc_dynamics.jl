struct CalciumDynamics{T<:Number} <: SpeciesDynamics{Calcium}
    τCa::T
    Ca∞::T
    calc_multiplier::T # f * area = 14.96 * 0.0628
end

function (p::CalciumDynamics)(n::Neuron, vars, varnames, currents)
    τCa, Ca∞, calc_multiplier, Cₘ = get_parameters(p)
    Ca = vars[findfirst(x -> x == Calcium, varnames)]
    D(Ca) ~ (1 / τCa) * (-Ca + Ca∞ + (calc_multiplier * currents / Cₘ))
end

get_parameters(::CalciumDynamics) = @parameters τCa, Ca∞, calc_multiplier, Cₘ
default_params(p::CalciumDynamics, n::Neuron, vars, varnames) = Dict(
    get_parameters(p) .=> (p.τCa, p.Ca∞, p.calc_multiplier, capacitance(n.geometry))
)

struct CaReversalDynamics <: SpeciesDynamics{Calcium} end

function (p::CaReversalDynamics)(n::Neuron, vars, varnames, currents)
    Ca = vars[findfirst(x -> x == Calcium, varnames)]
    ECa = vars[findfirst(x -> x == Reversal{Calcium}, varnames)]
    return ECa ~ (500.0) * (8.6174e-5) * (283.15) * (log(max((3000.0 / Ca), 0.001)))
end
