"""
Constructing channels

needs the following fields:
name
dynamics
defaults


Needs the following functions:

sensed()
actuated()
tagged_internal_variables()
untagged_internal_variables()

is_dynamic() which defaults to the set of keys for channel.dynamics

THE PROBLEM IS THE KCA CHANNEL

"""

function untagged_internal_variables(b::BasicChannel)
    iter = Iterators.flatten(zip(keys(b.dynamics), keys(b.defaults)))
    return filter(x -> typeof(x) <: UntrackedQuantity, unique(iter))
end


struct BasicSingleIonChannel{I<:Species, F<:Function, N, Q1 <: Quantity, Q2<:Quantity} <: BasicChannel
    name::Symbol
    ion::I
    dynamics::Dict{Q1,F}
    defaults::Dict{Q2,N}
end

sensed(b::BasicSingleIonChannel) = Set((Voltage(), Reversal{typeof(b.ion)}(), b.ion))
actuated(b::BasicSingleIonChannel) = Set((Current{Voltage}(), Current{typeof(b.ion)}()))
tagged_internal_variables(b::BasicSingleIonChannel) = Set((Conductance{typeof(b.ion)}(),Conductance{Voltage}() ))



struct BasicMultipleIonChannel{F<:Function,N,Q1<:Quantity,Q2<:Quantity,Q3<:Quantity,Q4<:Quantity} <: BasicChannel
    name::Symbol
    sensed::Set{Q1}
    actuated::Set{Q2}
    dynamics::Dict{Q3,F}
    defaults::Dict{Q4,N}
end


# sensed(b::BasicMultipleIonChannel) = b.sensed

sensed(b::BasicMultipleIonChannel) = Set((Voltage(), b.sensed...)) # NOT EQUAL TO B.SENSED
actuated(b::BasicMultipleIonChannel) = Set(((Current{typeof(el)}() for el in b.actuated)..., Current{Voltage}()))
tagged_internal_variables(b::BasicMultipleIonChannel) = Set(Conductance{typeof(i)}() for i in b.actuated if typeof(i)<:Species) 

# Set(map(i -> Conductance{typeof(i)}(), b.sensed)...,)




"""
KCa channel is calsium sensitive potassium


Sensed: reversal(actuated ion), Voltage, sensed ions.
"""