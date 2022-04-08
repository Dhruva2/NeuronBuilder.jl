abstract type Component end
abstract type IonChannel <: Component end
abstract type RegIonChannel <: Component end

struct RegIon{F<:AbstractFloat,I<:IonChannel} <: RegIonChannel
    ionch::I
    Rion::F
    τion::Symbol
end

abstract type Synapse <: Component end
struct EmptyConnection <: Synapse end
abstract type Compartment <: Component end

struct Soma{F<:AbstractFloat} <: Compartment
    initial_states::Dict{Symbol,F}
    parameters::Dict{Symbol,F}
end


function Soma_Requirements(channels::Vector{T}) where {T<:Component}
    mapreduce((x, y) -> unique!([x..., y...]), channels) do channel
        external_params(channel)
    end
end

"""
    Prinz_conversion(s::Soma)
This function converts bare parameters in units of milliSiemens/cm2 (as used by the prinz paper) to conform with units of microSiemens/mm2

"""
Prinz_conversion(s::Soma) = s.parameters[:Cₘ] / s.parameters[:area]

"""
    Liu_conversion(s::Soma)
This function converts bare parameters in units of microSiemens/nanoFarads (as used by the Liu paper) to conform with units of microSiemens/mm2

"""
Liu_conversion(s::Soma) = s.parameters[:Cₘ]

"""
    syn_conv_factor(s::Soma)
This function converts bare synaptic parameters in units of nanoSiemens*area (as used by the prinz paper, where area is baked in the synaptic conductances) to conform with units of microSiemens/mm2

"""
syn_conv_factor(s::Soma) = 1e-3 / s.parameters[:area]^2 #this gives μS/mm^2 

struct ComponentSystem{C<:Component}
    c::C
    sys::ODESystem
end


