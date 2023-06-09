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

sensed(b::BasicSingleIonChannel) = (Voltage(), Reversal{typeof(b.ion)}())
actuated(b::BasicSingleIonChannel) = (Current{Voltage}(), Current{typeof(b.ion)}())
tagged_internal_variables(b::BasicSingleIonChannel) = (Conductance{typeof(b.ion)}(),)



### delete this:
has_dynamics(::BasicSingleIonChannel, ::Voltage) = true
has_dynamics(::BasicSingleIonChannel, ::Calcium) = true
has_dynamics(::BasicSingleIonChannel, ::Reversal{S}) where S = false


struct BasicMultipleIonChannel{F<:Function,N,Q1<:Quantity,Q2<:Quantity,Q3<:Quantity,Q4<:Quantity} <: BasicChannel
    name::Symbol
    sensed::Vector{Q1}
    actuated::Vector{Q2}
    dynamics::Dict{Q3,F}
    defaults::Dict{Q4,N}
end


sensed(b::BasicMultipleIonChannel) = (Voltage(), b.sensed..., [Reversal{typeof(el)}() for el in b.sensed]...)
actuated(b::BasicMultipleIonChannel) = ([Current{typeof(el)}() for el in b.actuated]..., Current{Voltage}())
tagged_internal_variables(b::BasicMultipleIonChannel) = (map(i -> Conductance{typeof(i)}(), b.sensed)...,)




"""
KCa channel is calsium sensitive potassium

gatings:
m^4
m^3 h
m
no gate (leak)
"""