mutable struct GABA_A{T<:AbstractFloat} <: Synapse
    ḡGABA_A::T
    s::T
    Eₛ::T
    τ_R::T
    τ_D::T
end

GABA_A(x) = GABA_A(x, 0.0, -80.0, 0.5, 10.)
syn_current(::GABA_A, sys::ODESystem) = sys.IGABA_A

function channel_dynamics(channel::GABA_A, Vpre, Vpost)
    states = @variables s(t) IGABA_A(t)
    parameters = @parameters ḡGABA_A

    eqs = [
        D(s) ~ 0.5 * (1 + tanh(Vpre / 10.)) * ((1 - s) / channel.τ_R) - s / channel.τ_D,
        IGABA_A ~ -ḡGABA_A * s * (Vpost - channel.Eₛ)
        ]
    current = [eqs[2]]

    defaultmap = [s => channel.s, ḡGABA_A => channel.ḡGABA_A]

    return eqs, states, parameters, current, defaultmap
end

mutable struct AMPA{T<:AbstractFloat} <: Synapse
    ḡAMPA::T
    s::T
    Eₛ::T
    τ_R::T
    τ_D::T
end

AMPA(x) = AMPA(x, 0.0, 0.0, 0.2, 2.)
syn_current(::AMPA, sys::ODESystem) = sys.IAMPA

function channel_dynamics(channel::AMPA, Vpre, Vpost)
    states = @variables s(t) IAMPA(t)
    parameters = @parameters ḡAMPA

    eqs = [
        D(s) ~ 0.5 * (1 + tanh(Vpre / 10.)) * ((1 - s) / channel.τ_R) - s / channel.τ_D,
        IAMPA ~ -ḡAMPA * s * (Vpost - channel.Eₛ)
        ]
    current = [eqs[2]]

    defaultmap = [s => channel.s, ḡAMPA => channel.ḡAMPA]

    return eqs, states, parameters, current, defaultmap
end