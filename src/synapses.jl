
get_pre(::Synapse, sys::ODESystem) = sys.V
get_post(::Synapse, sys::ODESystem) = sys.V

my_pre(syn_sys::ComponentSystem) = syn_sys.sys.Vpre
my_post(syn_sys::ComponentSystem) = syn_sys.sys.Vpost
post_connector(::Synapse) = :Isyn

mutable struct Chol{T<:AbstractFloat} <: Synapse
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


mutable struct Glut{T<:AbstractFloat} <: Synapse
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


function channel_dynamics(ch::Chol, Vpre, Vpost)
    states = @variables s(t) IChol(t)
    parameters = @parameters ḡChol
    eqs = [D(s) ~ (1 / τs(ch, Vpre)) * (s̄(ch, Vpre) - s),
        IChol ~ -ḡChol * s * (Vpost - ch.Eₛ)]
    current = [eqs[2]]
    defaultmap = [s => ch.s, ḡChol => ch.ḡChol]
    return eqs, states, parameters, current, defaultmap
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
