### Basic components and subtypes ###

abstract type Species end
abstract type SpeciesProperty{S<:Species} end
abstract type SpeciesDynamics{F<:Species} end
abstract type Ion <: Species end
abstract type PseudoIon <: Ion end

"""
Not including voltage in tagging ion channels. Flows are all charged particles. So they will always include voltage.
"""

struct Voltage <: Species end
struct Sodium <: Ion end
struct Potassium <: Ion end
struct Calcium <: Ion end
struct Proton <: Ion end
struct Leak <: PseudoIon end



abstract type Reversal{I<:Ion} <: SpeciesProperty{I} end
abstract type Current{I<:Ion} <: SpeciesProperty{I} end
abstract type Conductance{I<:Ion} <: SpeciesProperty{I} end
abstract type mRNA{I<:Ion} <: SpeciesProperty{I} end

shorthand_name(::Type{Voltage}) = :V
shorthand_name(::Type{Sodium}) = :Na
shorthand_name(::Type{Potassium}) = :K
shorthand_name(::Type{Calcium}) = :Ca
shorthand_name(::Type{Proton}) = :H
shorthand_name(::Type{Leak}) = :Leak
# shorthand_name(::Type{MixedIon} )
shorthand_name(::Type{Reversal{T}}) where {T} = Symbol(:E, shorthand_name(T))
shorthand_name(::Type{Current{T}}) where {T} = Symbol(:I, shorthand_name(T))
shorthand_name(::Type{Conductance{T}}) where {T} = Symbol(:g, shorthand_name(T))
shorthand_name(::Type{mRNA{T}}) where {T} = Symbol(:mRNA_, shorthand_name(T))
shorthand_name(x::Type{Tuple{T,R}}) where {T,R} = shorthand_name.(x.types)