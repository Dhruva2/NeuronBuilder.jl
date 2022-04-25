#################### Transient sodium current ###############################

struct Na{F,D<:Real} <: FlowChannel(Sodium)
    g::D
    m::F
    h::F
end

Na(x) = Na(x, 0.0, 0.0)

αₘ(::Na, V) = 0.1(35.0 + V) / (-exp(-(35.0 + V) / 10.0) + 1.0)
βₘ(::Na, V) = 4exp(-(V + 60) / 18.0)
m∞(ch::Na, V) = αₘ(ch,V)/(αₘ(ch,V)+βₘ(ch,V))

αₕ(::Na, V) = 0.07exp(-(V + 60.0) / 20.0)
βₕ(::Na, V) = 1.0 / (1 + exp((-30.0 - V) / 10.0))
h∞(ch::Na, V) = αₕ(ch,V)/(αₕ(ch,V)+βₕ(ch,V))

function (ch::Na)(n::Neuron; name=get_name(ch))

    I, E, V = instantiate_hooks(n, ch, currents, reversals, voltage)
    g, = instantiate_parameters(ch, conductances)
    @variables m(t) h(t)

    states, params = vardivide(V, m, h, I, g, E)
    eqs = [
        D(m) ~ αₘ(ch,V) * (1 - m) - βₘ(ch,V) * m,
        D(h) ~ αₕ(ch,V) * (1 - h) - βₕ(ch,V) * h,
        I ~ g * m^3 * h * (E - V)
    ]

    V₀ = n.somatic_parameters[Voltage]
    defaultmap = [m => m∞(ch,V₀), h => h∞(ch,V₀), g => ch.g]

    return ODESystem(eqs, t, states, params; defaults=defaultmap, name=name)
end

#################### Delayed rectifier potassium current ######################

struct K{F,D<:Real} <: FlowChannel(Potassium)
    g::D
    m::F
end

K(x) = K(x, 0.0)

αₘ(::K, V) = 0.01(V + 50.0) / (1 - exp(-(50.0 + V) / 10.0))
βₘ(::K, V) = 0.125exp(-(V + 60) / 80.0)
m∞(ch::K, V) = αₘ(ch,V)/(αₘ(ch,V)+βₘ(ch,V))

function (ch::K)(n::Neuron; name=get_name(ch))
    I, E, V = instantiate_hooks(n, ch, currents, reversals, voltage)
    g, = instantiate_parameters(ch, conductances)
    @variables m(t)

    eqs = [
        D(m) ~ αₘ(ch,V) * (1 - m) - βₘ(ch,V) * m,
        I ~ g * m^4 * (E - V)
    ]

    states, params = vardivide(V, E, m, I, g)

    V₀ = n.somatic_parameters[Voltage]
    defaultmap = [m => m∞(ch,V₀), g => ch.g]

    return ODESystem(eqs, t, states, params; defaults=defaultmap, name=name)
end

#################### Leak current #########################

struct leak{D<:Real} <: FlowChannel(Leak)
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


