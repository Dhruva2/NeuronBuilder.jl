abstract type Component end

abstract type IonChannel <: Component end
abstract type RegIonChannel <: Component end

abstract type Synapse <: Component end

abstract type Compartment <: Component end

mutable struct Soma <: Compartment
    initial_states::Dict{Symbol,Float64}
    parameters::Dict{Symbol,Float64}
    hooks::Int64
end

struct ComponentSystem
    c::Component
    sys::ODESystem
end


