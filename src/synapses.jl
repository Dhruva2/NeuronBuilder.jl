"""
### What to do 

- Want to generalise synapses to allow for:
1. actuation of ionic currents, as well as voltage
2. sensing of general species, as well as voltage 

How? One proposal:

- Each synapse builds an ODESystem similar to channels as now.
- The synapse abstract supertype has {pre_sensed, post_sensed, actuated}: Nested tuples of the relevant species. 

- build_neuron instantiates current variables for each synapse, indexed (in their name) by:
    - synapse number 
    - actuated species

- A new function: 

function connecting_equations(synapse_system, pre_neuron_system, post_neuron_system)

end



- want to preserve synapses being at the same level in the system hierarchy as the neurons.

### Actual necessary connections:

neuron:
tracked fluxes += get_actuator(synapse, thing) for thing in tracked_species
eg I(t) += Chol.Imix 
I_na(t) += Chol.Ina

synapse:
for el in pre_sensed, voltage
    synapse.el_pre ~ neuron_pre.sensed
end

for el in post_sensed
    synapse.el_post ~ neuron_post.sensed
end


function get_connecting_equations(synapse_sys, neuron_sys)
    return pre_connections, post_connections
end

Neuron has to have a hook, depending on synapse number (as before). But hook is for each actuated speices from the synapse, if the dynamics of that species are tracked by the neuron.

"""

pre_sensed(::Synapse{Spre,Spost,A}) where {Spre,Spost,A} = typeflatten(Spre)
post_sensed(::Synapse{Spre,Spost,A}) where {Spre,Spost,A} = typeflatten(Spost)
actuated(::Synapse{Spre,Spost,A}) where {Spre,Spost,A} = typeflatten(A)

export pre_sensed, post_sensed

get_pre(::Synapse, sys::ODESystem) = sys.V
get_post(::Synapse, sys::ODESystem) = sys.V

my_pre(::Synapse, syn_sys::AbstractTimeDependentSystem) = syn_sys.Vpre
my_post(::Synapse, syn_sys::AbstractTimeDependentSystem) = syn_sys.Vpost
post_connector(::Synapse) = :Isyn

connecting_voltages(::Synapse, ::Neuron, ::Neuron) = @variables Vpre, Vpost


mutable struct Chol{T<:AbstractFloat} <: Synapse(Nothing, Nothing, MixedIon)
    ḡChol::T
    s::T
    Eₛ::T # mV
    k₋::T # ms
    Vth::T # mV
    δ::T   # mV
end

Chol(x) = Chol(x, 0.0, -80.0, 0.01, -35.0, 5.0)
s̄(syn::Chol, Vpre) = 1.0 / (1.0 + exp((syn.Vth - Vpre) / syn.δ))
τs(syn::Chol, Vpre) = (1.0 - s̄(syn, Vpre)) / syn.k₋
syn_current(::Chol, sys::ODESystem) = sys.IChol


mutable struct Glut{T<:AbstractFloat} <: Synapse(Nothing, Nothing, MixedIon)
    ḡGlut::T
    s::T
    Eₛ::T
    k₋::T
    Vth::T
    δ::T
end

Glut(x) = Glut(x, 0.0, -70.0, 0.025, -35.0, 5.0)
s̄(syn::Glut, Vpre) = 1.0 / (1.0 + exp((syn.Vth - Vpre) / syn.δ))
τs(syn::Glut, Vpre) = (1.0 - s̄(syn, Vpre)) / syn.k₋
syn_current(::Glut, sys::ODESystem) = sys.IGlut


function (ch::Chol)(n_pre::Neuron, n_post::Neuron; name=get_name(ch))
    I, = instantiate_variables(ch, currents)
    @variables s(t)
    Vpre, Vpost = connecting_voltages(ch, n_pre, n_post)
    @parameters ḡChol
    eqs = [D(s) ~ (1 / τs(ch, Vpre)) * (s̄(ch, Vpre) - s),
        I ~ -ḡChol * s * (Vpost - ch.Eₛ)]
    defaultmap = [s => ch.s, ḡChol => ch.ḡChol]
    return ODESystem(eqs, t, [I, s], parameters; defaults=defaultmap, name)
end

function channel_dynamics(ch::Glut, Vpre, Vpost)
    states = @variables s(t) IGlut(t)
    parameters = @parameters ḡGlut
    eqs = [D(s) ~ (1 / τs(ch, Vpre)) * (s̄(ch, Vpre) - s),
        IGlut ~ -ḡGlut * s * (Vpost - ch.Eₛ)]
    current = [eqs[2]]
    defaultmap = [s => ch.s, ḡGlut => ch.ḡGlut]
    return eqs, states, parameters, current, defaultmap
end

# function (syn::Synapse)(; name=get_name(syn))
#     @variables Vpre(t) Vpost(t)
#     eqs, states, parameters, current, defaultmap = channel_dynamics(syn, Vpre, Vpost)
#     sys = ODESystem(eqs, t, [Vpre, Vpost, states...], [parameters...]; observed=current,
#         name=name, defaults=defaultmap)
#     return sys
# end

# struct directed_synapse{S<:Synapse,N} <: Synapse
#     pre_n::N
#     post_n::N
#     syn::S
# end

struct TestSynapse2 <: Synapse(Calcium, Potassium, Chloride)
    g::Float64
end

export TestSynapse2

"""
1. get_sensor(;namespace=true) for each pre_sensed. and voltage. 
2. pre equations is easier: get_pre_sensed.(species) from the pre_neuron
3. 

end. return (ODESystem, pre_equations, post_equations)
"""
function (t::TestSynapse)(n_pre::Neuron, n_post::Neuron)

end