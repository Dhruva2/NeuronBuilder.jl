"""
regulated is a plasticity PlasticityRule


Plasticity rule can keep the original channel. 
Should be a component so it can be called as such from build_neuron
Does the needful for each regulated.
Stores extra plasticity-related variables separately

Callin a PlasticityRule() can then do the same as RegIonChannel, but using pr.channel, i.e. the channel stuff is separate from the plasticity rule in the namespace
"""


abstract type PlasticityRule end
component(::PlasticityRule) = PlasticityRule.component


struct OLearyCalciumReg{F<:AbstractFloat,I<:IonChannel} <: PlasticityRule
    component::I

end