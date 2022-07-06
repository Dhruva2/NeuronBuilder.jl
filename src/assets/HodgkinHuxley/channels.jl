αₙ(x) = 0.01(x + 50.0) / (1 - exp(
    -(50.0 + x) / 10.0
    )
)

βₙ(x) = 0.125exp(-(x + 60) / 80.0)

αₘ(x) = 0.1(35.0 + x) / (
    -exp(
        -(35.0 + x) / 10.0
    ) + 1.0
)

βₘ(x) = 4exp(
    -(x + 60) / 18.0
)

αₕ(x) = 0.07exp(-(x + 60.0) / 20.0);

βₕ(x) = 1.0 / (
    1 + exp(
        (-30.0 - x) / 10.0
    )
);

τₙ(x) = 1.0 / (αₙ(x) + βₙ(x));

"""
    Sensors: Reversal{Sodium}
    Actuators: Sodium
    Fields: Default values you might want
"""
struct Na{F,D} <: FlowChannel(Sodium)
    g::D
    m::F
    h::F
end

"""
    Default constructor
"""
Na(x) = Na(x, 0., 0.)

function (channel::Na)(neuron::Neuron; name = get_name(channel))

    # Instantiate any variables/parameters that need to be called programmatically
    I, E, V, g = instantiate_hooks(neuron, channel,
    currents, reversals, voltage, conductances)

    # Instantiate internal variables as you like
    @variables m(t) h(t)

    # Designate respective states / parameters
    states, params = vardivide(V, m, h, I, g, E)
    eqs = [
        D(m) ~ αₘ(V) * (1 - m) - βₘ(V) * m,
        D(h) ~ αₕ(V) * (1 - h) - βₕ(V) * h,
        I ~ g * m^3 * h * (E - V)
    ]

    defaultmap = [m => channel.m, h => channel.h, g => channel.g]

    return ODESystem(eqs, t, states, params; defaults = defaultmap, name = name)
end

struct K{F,D} <: FlowChannel(Potassium)
    g::D
    n::F
end

K(x) = K(x, 0.0)

function (ch::K)(neur::Neuron; name=get_name(ch))

    I, E, V = instantiate_hooks(neur, ch, currents, reversals, voltage)
    g, = instantiate_parameters(ch, conductances)
    @variables n(t)

    states, params = vardivide(V, n, I, g, E)
    eqs = [
        D(n) ~ αₙ(V) * (1 - n) - βₙ(V) * n,
        I ~ g * n^4 * (E - V)
    ]

    defaultmap = [n => ch.n, g => ch.g]
    return ODESystem(eqs, t, states, params; defaults=defaultmap, name=name)
end

struct Inp{F} <: FlowChannel{Tuple{Nothing},Tuple{Leak}}
    f::F
end

function (ch::Inp)(neur::Neuron; name=get_name(ch))
    I, V = instantiate_hooks(neur, ch, currents, voltage)

    eqs = [I ~ ch.f(t)]
    states, params = vardivide(V, I)
    return ODESystem(eqs, t, states, params; defaults=[I => ch.f(0)], name)
end

struct leak{D} <: FlowChannel(Leak)
    g::D
end

function (ch::leak)(n::Neuron; name=get_name(ch))
    I, E, V = instantiate_hooks(n, ch, currents, reversals, voltage)
    g, = instantiate_parameters(ch, conductances)

    eqs = [I ~ g * (E - V)]

    states, params = vardivide(V, E, I, g)
    defaultmap = [g => ch.g]

    return ODESystem(eqs, t, states, params; defaults=defaultmap, name=name)
end

struct Input{F} <: FlowChannel{Tuple{Nothing},Tuple{Leak}}
    f::F
end

function (ch::Input)(neur::Neuron; name=get_name(ch))
    I, V = instantiate_hooks(neur, ch, currents, voltage)

    eqs = [I ~ ch.f(t)]
    states, params = vardivide(V, I)
    return ODESystem(eqs, t, states, params; defaults=[I => ch.f(0)], name)
end