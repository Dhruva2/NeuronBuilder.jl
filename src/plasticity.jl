"""
regulated is a plasticity PlasticityRule


Plasticity rule can keep the original channel. 
Should be a component so it can be called as such from build_neuron
Does the needful for each regulated.
Stores extra plasticity-related variables separately

Callin a PlasticityRule() can then do the same as RegIonChannel, but using pr.channel, i.e. the channel stuff is separate from the plasticity rule in the namespace
"""


component(::PlasticityRule) = PlasticityRule.component


struct IntegralCalciumReg{F<:AbstractFloat,I<:IonChannel} <: PlasticityRule
    component::I
    Rion::F
    τion::Symbol
end

"""
1. get name of parameters to turn into states
2. create new states out of them
3. append equations
1. and 2. might be independent of plasticity rule
"""
function add_plasticity_dynamics(I::IntegralCalciumReg)
    @variables V(t) Ca(t)
    eqs, states, parameters, current, defaultmap = channel_dynamics(I.component, V, Ca)
    
end


function (pr::PlasticityRule)()

end

