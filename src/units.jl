#unit handling
using Unitful
using Unitful.DefaultSymbols
import Unitful: mV, mS, cm, mm, nA, mA, µA, ms, nF, μM
import Unitful: Time, Molarity, ElectricalConductance, Area
import Unitful: Voltage as voltage
import Unitful: Current as current

#NOTE Iapp is *specific* current in nA / mm^2

@derived_dimension SpecificCapacitance 𝐈^2 * 𝐋^-4 * 𝐌^-1 * 𝐓^4 # capacitance per unit area. SI base units
function Cm(Cm::Unitful.Quantity=10.0nF / mm^2)
    !(typeof(Cm) <: SpecificCapacitance) && throw("Specific capacitance must have units of Farads / squared length")
    uconvert(nF / mm^2, Cm)
end

function Area(area::Unitful.Quantity=0.0628mm^2)
    (dimension(area) != 𝐋^2) && throw("Area must have units of squared length")
    uconvert(mm^2, area)
end

function check_units(dic::Dict{Symbol,Q}) where {Q<:Unitful.Quantity}
    haskey(dic, :Ca∞) && (get(dic, :Ca∞, 0) |> q -> typeof(q) <: Unitful.Molarity ? dic[:Ca∞] = ustrip(Real, μM, q) : throw("Check calcium units are μM"))
    haskey(dic, :area) && (get(dic, :area, 0) |> q -> typeof(q) <: Unitful.Area ? dic[:area] = ustrip(Real, mm^2, q) : throw("Check area units are mm^2"))
    haskey(dic, :Cₘ) && (get(dic, :Cₘ, 0) |> q -> typeof(q) <: SpecificCapacitance ? dic[:Cₘ] = ustrip(Real, nF / mm^2, q) : throw("Check specific capacitance units are nF/mm^2"))
    haskey(dic, :τCa) && (get(dic, :τCa, 0) |> q -> typeof(q) <: Unitful.Time ? dic[:τCa] = ustrip(Real, ms, q) : throw("Check time units are ms"))
    haskey(dic, :Iapp) && (get(dic, :Iapp, 0) |> q -> typeof(q) <: current ? dic[:Iapp] = ustrip(Real, nA / mm^2, q) : throw("Check current units are nA/mm^2"))

    return dic
    #returns bare values in the right units
end

function check_units(dic::Dict{DataType,Q}) where {Q<:Unitful.Quantity}
    haskey(dic, Voltage) && (get(dic, Voltage, 0) |> q -> typeof(q) <: voltage ? dic[Voltage] = ustrip(Real, mV, q) : throw("Check voltage units are mV"))
    haskey(dic, Calcium) && (get(dic, Calcium, 0) |> q -> typeof(q) <: Unitful.Molarity ? dic[Calcium] = ustrip(Real, μM, q) : throw("Check calcium units are μM"))
    return dic
end

ICS(dict) = check_units(dict)
Params(dict) = check_units(dict)

function Reversals(dic::Dict{DataType,Q}) where {Q<:Unitful.Quantity}
    Dict((kv[1], ustrip(Real, mV, kv[2])) for kv in dic)
end