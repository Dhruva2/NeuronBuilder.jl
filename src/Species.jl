

abstract type Quantity end
abstract type TrackedQuantity <: Quantity end
abstract type Species <: TrackedQuantity end
abstract type SpeciesProperty{S<:Species} <: TrackedQuantity end
abstract type SpeciesDynamics{F<:Species} <:TrackedQuantity end
abstract type Ion <: Species end
abstract type PseudoIon <: Ion end

struct UntrackedQuantity <: Quantity
name::Symbol
end

"""
Not including voltage in tagging ion channels. Flows are all charged particles. So they will always include voltage.
"""

struct Voltage <: Species end
struct Sodium <: Ion end
struct Potassium <: Ion end
struct Calcium <: Ion end
struct Proton <: Ion end
struct Leak <: PseudoIon end



struct Reversal{S<:Species} <: SpeciesProperty{S} end
struct Current{S<:Species} <: SpeciesProperty{S} end
struct Conductance{S<:Species} <: SpeciesProperty{S} end
struct mRNA{I<:Ion} <: SpeciesProperty{I} end
struct Capacitance{S<:Species} <: SpeciesProperty{S} end

### all of these shorthand/prefix names should be unique!!
shorthand_name(::Voltage) = :V
shorthand_name(::Sodium) = :Na
shorthand_name(::Potassium) = :K
shorthand_name(::Calcium) = :Ca
shorthand_name(::Proton) = :H
shorthand_name(::Leak) = :Leak


prefix_name(::Reversal) = :E 
prefix_name(::Current) = :I 
prefix_name(::Conductance) = :g
prefix_name(::mRNA) = :mRNA
prefix_name(::Capacitance) = :C


shorthand_name(S::SpeciesProperty{I}) where {I} = Symbol(prefix_name(S), shorthand_name(I()))
shorthand_name(S::SpeciesProperty{Voltage}) = prefix_name(S) 

shorthand_name(u::UntrackedQuantity) = u.name


# abstract type QuantityType end 
# struct Flux <: QuantityType end
# struct Fluence <: QuantityType end

# variable_type(::Species) = Fluence()
# variable_type(::Reversal) = Fluence() # doesn't determine dynamics. eg reversal an be dynamic
# variable_type(::Conductance) = Fluence()
# variable_type(::mRNA) = Fluence()
# variable_type(::Current) = Flux()

