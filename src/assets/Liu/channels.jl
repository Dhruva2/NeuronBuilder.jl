#################### NaV ###############################

struct Na{F,D<:Real} <: FlowChannel(Sodium)
    gNa::D
    mNa::F
    hNa::F
end

Na(x) = Na(x, 0.0, 0.0)
m∞(::Na, V) = 1.0 / (1.0 + exp((V + 25.5) / -5.29))
h∞(::Na, V) = 1.0 / (1.0 + exp((V + 48.9) / 5.18))
τm(::Na, V) = 1.32 - 1.26 / (1 + exp((V + 120.0) / -25.0))
τh(::Na, V) = (0.67 / (1.0 + exp((V + 62.9) / -10.0))) * (1.5 + 1.0 / (1.0 + exp((V + 34.9) / 3.6)))

function (ch::Na)(n::Neuron;name=get_name(ch))

    I, E, V = instantiate_hooks(n, ch, currents, reversals, voltage)
    g, = instantiate_parameters(ch, conductances)
    @variables mNa(t) hNa(t)

    states, params = vardivide(V, mNa, hNa, I, g, E)
    eqs = [D(mNa) ~ (1 / τm(ch, V)) * (m∞(ch, V) - mNa),
        D(hNa) ~ (1 / τh(ch, V)) * (h∞(ch, V) - hNa),
        I ~ g * mNa^3 * hNa * (E - V)]

    defaultmap = [mNa => m∞(ch, n.somatic_parameters[Voltage]), hNa => h∞(ch, n.somatic_parameters[Voltage]), g => ch.gNa]
    return ComponentSystem(ch, ODESystem(eqs, t, states, params; defaults=defaultmap, name=name))
end

#################### Slow calcium current #############################
#doesn't have to sense calcium. only ECa needs ca sensed but  
struct CaS{F,D<:Real} <: FlowChannel(Calcium)
    gCaS::D
    mCaS::F
    hCaS::F
end

CaS(x) = CaS(x, 20.0, 0.0)
m∞(::CaS, V) = 1.0 / (1.0 + exp((V + 33.0) / -8.1))
h∞(::CaS, V) = 1.0 / (1.0 + exp((V + 60.0) / 6.2))
τm(::CaS, V) = 1.4 + 7.0 / (exp((V + 27.0) / 10.0) + exp((V + 70.0) / -13.0))
τh(::CaS, V) = 60.0 + 150.0 / (exp((V + 55.0) / 9.0) + exp((V + 65.0) / -16.0))

function (ch::CaS)(n::Neuron;name=get_name(ch))
    I, E, V = instantiate_hooks(n, ch, currents, reversals, voltage)
    g, = instantiate_parameters(ch, conductances)
    @variables mCaS(t) hCaS(t) V(t)

    eqs = [D(mCaS) ~ (1 / τm(ch, V)) * (m∞(ch, V) - mCaS),
        D(hCaS) ~ (1 / τh(ch, V)) * (h∞(ch, V) - hCaS),
        I ~ g * mCaS^3 * hCaS * (E - V)]

    states, params = vardivide(V, E, mCaS, hCaS, I, g)
    defaultmap = [mCaS => m∞(ch, n.somatic_parameters[Voltage]), hCaS => h∞(ch, n.somatic_parameters[Voltage]), g => ch.gCaS]

    return ComponentSystem(ch, ODESystem(eqs, t, states, params; defaults=defaultmap, name=name))
end

#################### Transient calcium current ######################

struct CaT{F,D<:Real} <: FlowChannel(Calcium)
    gCaT::D
    mCaT::F
    hCaT::F
end

CaT(x) = CaT(x, 0.0, 0.0)
m∞(::CaT, V) = 1.0 / (1.0 + exp((V + 27.1) / -7.2))
h∞(::CaT, V) = 1.0 / (1.0 + exp((V + 32.1) / 5.5))
τm(::CaT, V) = 21.7 - 21.3 / (1.0 + exp((V + 68.1) / -20.5));
τh(::CaT, V) = 105.0 - 89.8 / (1.0 + exp((V + 55.0) / -16.9));

function (ch::CaT)(n::Neuron;name=get_name(ch))
    I, E, V = instantiate_hooks(n, ch, currents, reversals, voltage)
    g, = instantiate_parameters(ch, conductances)

    @variables mCaT(t) hCaT(t)
    eqs = [D(mCaT) ~ (1 / τm(ch, V)) * (m∞(ch, V) - mCaT),
        D(hCaT) ~ (1 / τh(ch, V)) * (h∞(ch, V) - hCaT),
        I ~ g * mCaT^3 * hCaT * (E - V)]

    states, params = vardivide(V, E, mCaT, hCaT, I, g)
    defaultmap = [mCaT => m∞(ch, n.somatic_parameters[Voltage]), hCaT => h∞(ch, n.somatic_parameters[Voltage]), g => ch.gCaT]

    return ComponentSystem(ch, ODESystem(eqs, t, states, params; defaults=defaultmap, name=name))
end

#################### A-type potassium current #########################

struct Ka{F,D<:Real} <: FlowChannel(Potassium)
    gKa::D
    mKa::F
    hKa::F
end

Ka(x) = Ka(x, 0.0, 0.0)
m∞(::Ka, V) = 1.0 / (1.0 + exp((V + 27.2) / -8.7))
h∞(::Ka, V) = 1.0 / (1.0 + exp((V + 56.9) / 4.9))
τm(::Ka, V) = 11.6 - 10.4 / (1.0 + exp((V + 32.9) / -15.2))
τh(::Ka, V) = 38.6 - 29.2 / (1.0 + exp((V + 38.9) / -26.5))

function (ch::Ka)(n::Neuron;name=get_name(ch))
    I, E, V = instantiate_hooks(n, ch, currents, reversals, voltage)
    g, = instantiate_parameters(ch, conductances)

    @variables mKa(t) hKa(t)
    eqs = [D(mKa) ~ (1 / τm(ch, V)) * (m∞(ch, V) - mKa),
        D(hKa) ~ (1 / τh(ch, V)) * (h∞(ch, V) - hKa),
        I ~ g * mKa^3 * hKa * (E - V)]

    states, params = vardivide(V, E, mKa, hKa, I, g)
    defaultmap = [mKa => m∞(ch, n.somatic_parameters[Voltage]), hKa => h∞(ch, n.somatic_parameters[Voltage]), g => ch.gKa]

    return ComponentSystem(ch, ODESystem(eqs, t, states, params; defaults=defaultmap, name=name))
end

################### Calcium-activated potassium current ########

struct KCa{F,D<:Real} <: FlowChannel(Calcium, Potassium)
    gKCa::D
    mKCa::F
end

KCa(x) = KCa(x, 0.0)
m∞(::KCa, V, Ca) = (Ca / (Ca + 3.0)) / (1.0 + exp((V + 28.3) / -12.6));
τm(::KCa, V) = 90.3 - 75.1 / (1.0 + exp((V + 46.0) / -22.7));

function (ch::KCa)(n::Neuron;name=get_name(ch))
    I, Ca, E, V = instantiate_hooks(n, ch, currents, sensed_ions, reversals, voltage)
    g, = instantiate_parameters(ch, conductances)

    @variables mKCa(t)
    eqs = [D(mKCa) ~ (1 / τm(ch, V)) * (m∞(ch, V, Ca) - mKCa),
        I ~ (g * mKCa^4) * (E - V)]

    states, params = vardivide(V, Ca, E, mKCa, I, g)
    defaultmap = [mKCa => m∞(ch, n.somatic_parameters[Voltage], n.somatic_parameters[Calcium]), g => ch.gKCa]

    return ComponentSystem(ch, ODESystem(eqs, t, states, params; defaults=defaultmap, name=name))
end

#################### Delayed rectifier potassium current ######################
struct Kdr{F,D<:Real} <: FlowChannel(Potassium)
    gKdr::D
    mKdr::F
end

Kdr(x) = Kdr(x, 0.0)
m∞(::Kdr, V) = 1.0 / (1.0 + exp((V + 12.3) / -11.8));
τm(::Kdr, V) = 7.2 - 6.4 / (1.0 + exp((V + 28.3) / -19.2));

function (ch::Kdr)(n::Neuron;name=get_name(ch))
    I, E, V = instantiate_hooks(n, ch, currents, reversals, voltage)
    g, = instantiate_parameters(ch, conductances)
    @variables mKdr(t)

    eqs = [D(mKdr) ~ (1 / τm(ch, V)) * (m∞(ch, V) - mKdr),
        I ~ (g * mKdr^4) * (E - V)]

    states, params = vardivide(V, E, mKdr, I, g)
    defaultmap = [mKdr => m∞(ch, n.somatic_parameters[Voltage]), g => ch.gKdr]

    return ComponentSystem(ch, ODESystem(eqs, t, states, params; defaults=defaultmap, name=name))
end

#################### H current ####################

struct H{F,D<:Real} <: FlowChannel(Proton)
    gH::D
    mH::F
end

H(x) = H(x, 0.0)
m∞(::H, V) = 1.0 / (1.0 + exp((V + 70.0) / 6.0))
τm(::H, V) = (272.0 + 1499.0 / (1.0 + exp((V + 42.2) / -8.73)))

function (ch::H)(n::Neuron;name=get_name(ch))
    I, E, V = instantiate_hooks(n, ch, currents, reversals, voltage)
    g, = instantiate_parameters(ch, conductances)
    @variables mH(t)

    eqs = [D(mH) ~ (1 / τm(ch, V)) * (m∞(ch, V) - mH),
        I ~ g * mH * (E - V)]

    states, params = vardivide(V, E, mH, I, g)
    defaultmap = [mH => m∞(ch, n.somatic_parameters[Voltage]), g => ch.gH]

    return ComponentSystem(ch, ODESystem(eqs, t, states, params; defaults=defaultmap, name=name))
end

#################### Leak current #########################

struct leak{D<:Real} <: FlowChannel(Leak)
    gLeak::D
end

function (ch::leak)(n::Neuron;name=get_name(ch))
    I, E, V = instantiate_hooks(n, ch, currents, reversals, voltage)
    g, = instantiate_parameters(ch, conductances)

    eqs = [I ~ g * (E - V)]

    states, params = vardivide(V, E, I, g)
    defaultmap = [g => ch.gLeak]

    return ComponentSystem(ch, ODESystem(eqs, t, states, params; defaults=defaultmap, name=name))
end


