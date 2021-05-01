abstract type Conductance end
abstract type IonChannel <: Conductance end
abstract type Synapse <: Conductance end

"""
return unparameterised type as symbol
"""
function get_name(ch::Conductance) 
    Base.typename(ch |> typeof).name |> Symbol
end

"""
fallback option for channels without a calcium current
"""
calcium_current(ch::IonChannel, sys::ODESystem) = Num(0)
voltage_hook(V, ch::IonChannel, sys::ODESystem) = V ~ sys.V 
calcium_hook(Ca, ch::IonChannel, sys::ODESystem) = Ca ~ sys.Ca

"""
does the channel have a reversal
"""
get_reversal(c::Conductance) = nothing
#################### Channels ###############################