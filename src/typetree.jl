### Basic components and subtypes
abstract type Flower end
abstract type Ion <: Flower end
abstract type AbstractIon <: Ion end

"""
Not including voltage in tagging ion channels. Any flow channel. Flows are all charged particles. So they will always include voltage. Not going to consider flows of e.g. proteins
"""
struct Voltage <: Flower end
struct Sodium <: Ion end
struct Potassium <: Ion end
struct Calcium <: Ion end
struct Proton <: Ion end
struct Leak <: AbstractIon end

ionic(x) = x <: Ion
export ionic

shorthand_name(::Type{Voltage}) = :V
shorthand_name(::Type{Sodium}) = :Na
shorthand_name(::Type{Potassium}) = :K
shorthand_name(::Type{Calcium}) = :Ca
shorthand_name(::Type{Proton}) = :H
shorthand_name(::Type{Leak}) = :Leak
export shorthand_name


abstract type Component end
abstract type Compartment <: Component end
abstract type FlowChannel{Sensors,Actuators} <: Component end

# is_geometric(::Compartment{T}) = T
# is_geometric(::Component) = false

FlowChannel(T) = FlowChannel{Tuple{},T}
FlowChannel(S, A) = FlowChannel{S,A}

sensed(::FlowChannel{S,A}) where {S<:Tuple} where {A} = fieldtypes(S)
actuated(::FlowChannel{S,A}) where {S} where {A<:Tuple} = fieldtypes(A)

sensed(::FlowChannel{S,A}) where {S<:Flower} where {A} = (S,)
actuated(::FlowChannel{S,A}) where {S} where {A<:Flower} = (A,)

export sensed, actuated

flows(i::FlowChannel) = (sensed(i)..., actuated(i)...) |> unique
ions(i::FlowChannel) = filter(x -> x <: Ion, flows(i)) |> unique

struct TestChannel1 <: FlowChannel{Sodium,Potassium} end
struct TestChannel2 <: FlowChannel{Tuple{Sodium,Calcium},Potassium} end
struct TestChannel3 <: FlowChannel{Tuple{Sodium,Calcium},Tuple{Potassium,Proton}} end
struct TestChannel4 <: FlowChannel{Tuple{Voltage,Sodium},Tuple{Voltage,Potassium}} end
export TestChannel1, TestChannel2, TestChannel3, TestChannel4

sensedvars(i::FlowChannel) =
    map(sensed(i)) do thing
        Symbol(thing |> shorthand_name)
    end

# function reversals(i::FlowChannel)
    
# end

reversals(i::FlowChannel) =
     map(actuated(i)) do thing
        Symbol(:E, shorthand_name(thing))
    end 

currents(i::FlowChannel) =
    map(actuated(i)) do thing
        Symbol(:I, shorthand_name(thing))
    end


conductances(i::FlowChannel) =
    map(actuated(i)) do thing
        Symbol(:g, shorthand_name(thing))
    end


instantiate_variables(c::Component, args...) =
    map(args) do el
        map(el(c)) do f
            @variables $f(t)
        end |> Iterators.flatten |> collect
    end |> Iterators.flatten |> collect

instantiate_parameters(c::Component, args...) =
    map(args) do el
        map(el(c)) do f
            @parameters $f
        end |> Iterators.flatten |> collect
    end |> Iterators.flatten |> collect


instantiate_variables(v::Vector{Symbol}) = 
    map(v) do f
        @variables $f(t)
    end |> Iterators.flatten |> collect  




export reversals, currents, flows, conductances, instantiate_variables, instantiate_parameters, sensedvars


abstract type Synapse <: Component end
struct EmptyConnection <: Synapse end

### delete these asap
abstract type IonChannel <: Component end
abstract type CalciumChannel <: IonChannel end
### 



abstract type RegIonChannel <: Component end # delete eventually
abstract type PlasticityRule <: Component end

struct Soma{F<:AbstractFloat} <: Compartment
    initial_states::Dict{Symbol,F}
    parameters::Dict{Symbol,F}
end



struct ComponentSystem{C<:Component,S<:AbstractTimeDependentSystem}
    c::C
    sys::S
end

# flows(c::ComponentSystem) = flows(c.c)
# sensedvars(c::ComponentSystem) = sensedvars(c.c)
# reversals(c::ComponentSystem) = reversals(c.c)
# currents(c::ComponentSystem) = currents(c.c)
# conductances(c::ComponentSystem) = conductances(c.c)

# export getva


"""
delete when PlasticityRule overrides
"""
struct RegIon{F<:AbstractFloat,I<:IonChannel} <: RegIonChannel
    ionch::I
    Rion::F
    τion::Symbol
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



######### Geometries ##############

abstract type Geometry end
struct NoGeometry{C} <: Geometry
    capacitance::C
end

abstract type Neuron <: Compartment end
