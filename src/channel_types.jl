abstract type Component end

abstract type IonChannel <: Component end
abstract type RegIonChannel <: Component end

struct RegIon{F<:AbstractFloat,I<:IonChannel} <: RegIonChannel
    ionch::I
    Rion::F
    τion::Symbol
end

abstract type Synapse <: Component end

abstract type Compartment <: Component end

mutable struct Soma{F<:AbstractFloat} <: Compartment
    initial_states::Dict{Symbol,F}
    parameters::Dict{Symbol,F}
    hooks::Int64
end

Soma(init_states, params) = Soma(init_states, params, 0)
"""
    Prinz_conversion(s::Soma)
This function does blah
"""
Prinz_conversion(s::Soma) = s.parameters[:Cₘ] / s.parameters[:area]
Liu_conversion(s::Soma) = s.parameters[:Cₘ]
syn_conv_factor(s::Soma) = 1e-3 / s.parameters[:area]^2 #this gives μS/mm^2 

struct ComponentSystem{C<:Component}
    c::C
    sys::ODESystem
end


