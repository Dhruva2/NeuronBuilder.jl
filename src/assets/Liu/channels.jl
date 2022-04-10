const f = 14.96 # uM*nF/nA 
"""
eg IonicConductance(::Na) = :gNa
"""

# NeuronBuilder.IonicCurrents(i::IonChannel) = (Symbol(:I, nameof(i |> typeof)),)


# # CalciumFlows(c::CalciumChannel) = (Symbol(:I, nameof(c |> typeof), :_Ca),)

# NeuronBuilder.Flows(i::IonChannel) = (:V,)
# # NeuronBuilder.Flows(c::CalciumChannel) = (:V, :Ca)
# NeuronBuilder.Reversals(i::IonChannel) = (Symbol(:E, nameof(i |> typeof)),)

# flows = NeuronBuilder.Flows
# ionic_conductances = NeuronBuilder.IonicConductances
# ionic_currents = NeuronBuilder.IonicCurrents
# reversals = NeuronBuilder.Reversals


#################### NaV ###############################

struct Na{F,D<:Real} <: FlowChannel(Sodium)
    gNa::D
    mNa::F
    hNa::F
end
Na(x) = Na(x, 0., 0.)
m∞(::Na, V) = 1.0 / (1.0 + exp((V + 25.5) / -5.29))
h∞(::Na, V) = 1.0 / (1.0 + exp((V + 48.9) / 5.18))
τm(::Na, V) = 1.32 - 1.26 / (1 + exp((V + 120.0) / -25.0))
τh(::Na, V) = (0.67 / (1.0 + exp((V + 62.9) / -10.0))) * (1.5 + 1.0 / (1.0 + exp((V + 34.9) / 3.6)))

function (ch::Na)()
    @variables mNa(t) hNa(t) V(t)
    (I,) = instantiate_variables(ch, currents)
    (E,), (g,) = instantiate_parameters(ch, reversals, conductances)

    eqs = [D(mNa) ~ (1 / τm(ch, V)) * (m∞(ch, V) - mNa),
        D(hNa) ~ (1 / τh(ch, V)) * (h∞(ch, V) - hNa),
        I ~ g * mNa^3 * hNa * (E - V)]
    current = [eqs[3]]
    defaultmap = [mNa => ch.mNa, hNa => ch.hNa, g => ch.gNa]
    return ComponentSystem(ch, ODESystem(eqs, t, [V, mNa, hNa, I], [g, E]; observed=current, defaults=defaultmap, name=get_name(ch)))
end


#################### Slow calcium current #############################

struct CaS{F,D<:Real} <: FlowChannel(Calcium, Calcium)
    gCaS::D
    mCaS::F
    hCaS::F
end

CaS(x) = CaS(x, 20.0, 0.0)
ECa(::CaS, Ca) = (500.0) * (8.6174e-5) * (283.15) * (log(max((3000.0 / Ca), 0.001))) #same as conductor 
ExternalHooks(::CaS) = (:τCaS,)
m∞(::CaS, V) = 1.0 / (1.0 + exp((V + 33.0) / -8.1))
h∞(::CaS, V) = 1.0 / (1.0 + exp((V + 60.0) / 6.2))
τm(::CaS, V) = 1.4 + 7.0 / (exp((V + 27.0) / 10.0) + exp((V + 70.0) / -13.0))
τh(::CaS, V) = 60.0 + 150.0 / (exp((V + 55.0) / 9.0) + exp((V + 65.0) / -16.0))

function (ch::CaS)()
    @variables mCaS(t) hCaS(t) V(t)
    (I,), (Ca,), (E,) = instantiate_variables(ch, currents, sensedvars, reversals)
    (g,) = instantiate_parameters(ch, conductances)

    eqs = [D(mCaS) ~ (1 / τm(ch, V)) * (m∞(ch, V) - mCaS),
        D(hCaS) ~ (1 / τh(ch, V)) * (h∞(ch, V) - hCaS),
        E ~ ECa(ch, Ca),
        I ~ g * mCaS^3 * hCaS * (E - V)]
    current = [eqs[4]]
    defaultmap = [mCaS => ch.mCaS, hCaS => ch.hCaS, g => ch.gCaS]
    return ComponentSystem(ch, ODESystem(eqs, t, [V, Ca, E, mCaS, hCaS, I], [g]; observed=current, defaults=defaultmap, name=get_name(ch)))
end

#################### Transient calcium current ######################

struct CaT{F,D<:Real} <: FlowChannel(Calcium, Calcium)
    gCaT::D
    mCaT::F
    hCaT::F
end

CaT(x) = CaT(x, 0.0, 0.0)
ECa(::CaT, Ca) = (500.0) * (8.6174e-5) * (283.15) * (log(max((3000.0 / Ca), 0.001)))
ExternalHooks(::CaT) = (:τCaT,)
m∞(::CaT, V) = 1.0 / (1.0 + exp((V + 27.1) / -7.2))
h∞(::CaT, V) = 1.0 / (1.0 + exp((V + 32.1) / 5.5))
τm(::CaT, V) = 21.7 - 21.3 / (1.0 + exp((V + 68.1) / -20.5));
τh(::CaT, V) = 105.0 - 89.8 / (1.0 + exp((V + 55.0) / -16.9));

function (ch::CaT)()
    @variables mCaT(t) hCaT(t) V(t)
    (I,), (Ca,), (E,) = instantiate_variables(ch, currents, sensedvars, reversals)
    (g,) = instantiate_parameters(ch, conductances)

    eqs = [D(mCaT) ~ (1 / τm(ch, V)) * (m∞(ch, V) - mCaT),
        D(hCaT) ~ (1 / τh(ch, V)) * (h∞(ch, V) - hCaT),
        E ~ ECa(ch, Ca),
        I ~ g * mCaT^3 * hCaT * (E - V)]
    current = [eqs[4]]
    defaultmap = [mCaT => ch.mCaT, hCaT => ch.hCaT, g => ch.gCaT]
    return ComponentSystem(ch, ODESystem(eqs, t, [V, Ca, E, mCaT, hCaT, I], [g]; observed=current, defaults=defaultmap, name=get_name(ch)))
end

#################### A-type potassium current #########################

struct Ka{F,D<:Real} <: FlowChannel(Potassium)
    gKa::D
    mKa::F
    hKa::F
end

Ka(x) = Ka(x, 0.0, 0.0)
ExternalHooks(::Ka) = (:EK, :τKa)
m∞(::Ka, V) = 1.0 / (1.0 + exp((V + 27.2) / -8.7))
h∞(::Ka, V) = 1.0 / (1.0 + exp((V + 56.9) / 4.9))
τm(::Ka, V) = 11.6 - 10.4 / (1.0 + exp((V + 32.9) / -15.2))
τh(::Ka, V) = 38.6 - 29.2 / (1.0 + exp((V + 38.9) / -26.5))

function (ch::Ka)()
    (I,) = instantiate_variables(ch, currents)
    (E,), (g,) = instantiate_parameters(ch, reversals, conductances)
    states = @variables mKa(t) hKa(t) V(t)
    eqs = [D(mKa) ~ (1 / τm(ch, V)) * (m∞(ch, V) - mKa),
        D(hKa) ~ (1 / τh(ch, V)) * (h∞(ch, V) - hKa),
        I ~ g * mKa^3 * hKa * (E - V)]
    current = [eqs[3]]
    defaultmap = [mKa => ch.mKa, hKa => ch.hKa, g => ch.gKa]
    return ComponentSystem(ch, ODESystem(eqs, t, [V, mKa, hKa, I], [g, E]; observed=current, defaults=defaultmap, name=get_name(ch)))
end

################### Calcium-activated potassium current ########

struct KCa{F,D<:Real} <: FlowChannel(Calcium, Potassium)
    gKCa::D
    mKCa::F
end

KCa(x) = KCa(x, 0.0)
ExternalHooks(::KCa) = (:EK, :τKCa)
m∞(::KCa, V, Ca) = (Ca / (Ca + 3.0)) / (1.0 + exp((V + 28.3) / -12.6));
τm(::KCa, V) = 90.3 - 75.1 / (1.0 + exp((V + 46.0) / -22.7));

function (ch::KCa)()
    @variables mKCa(t) V(t)
    I, Ca = instantiate_variables(ch, currents, sensedvars)
    g, E = instantiate_parameters(ch, conductances, reversals)
    eqs = [
        D(mKCa) ~ (1 / τm(ch, V)) * (m∞(ch, V, Ca) - mKCa),
        I ~ (g * mKCa^4) * (E - V)
        ]
    current = [eqs[2]]
    defaultmap = [mKCa => ch.mKCa, g => ch.gKCa]
    return ComponentSystem(ch, ODESystem(eqs, t, [V, Ca, mKCa, I], [g, E]; observed=current, defaults=defaultmap, name=get_name(ch)))
end

#################### Delayed rectifier potassium current ######################
struct Kdr{F,D<:Real} <: FlowChannel(Potassium)
    gKdr::D
    mKdr::F
end

Kdr(x) = Kdr(x, 0.0)
m∞(::Kdr, V) = 1.0 / (1.0 + exp((V + 12.3) / -11.8));
τm(::Kdr, V) = 7.2 - 6.4 / (1.0 + exp((V + 28.3) / -19.2));

function (ch::Kdr)()
    states = @variables mKdr(t) V(t)
    (I,) = instantiate_variables(ch, currents)
    (E,), (g,) = instantiate_parameters(ch, reversals, conductances)
    eqs = [D(mKdr) ~ (1 / τm(ch, V)) * (m∞(ch, V) - mKdr),
        I ~ (g * mKdr^4) * (E - V)]
    current = [eqs[2]]
    defaultmap = [mKdr => ch.mKdr, g => ch.gKdr]
    return ComponentSystem(ch, ODESystem(eqs, t, [V, mKdr, I], [g, E]; observed=current, defaults=defaultmap, name=get_name(ch)))
end

#################### H current ####################

struct H{F,D<:Real} <: FlowChannel(Proton)
    gH::D
    mH::F
end

H(x) = H(x, 0.0)
m∞(::H, V) = 1.0 / (1.0 + exp((V + 70.0) / 6.0))
τm(::H, V) = (272.0 + 1499.0 / (1.0 + exp((V + 42.2) / -8.73)))

function (ch::H)()
    @variables mH(t) V(t)
    (I,) = instantiate_variables(ch, currents)
    (E,), (g,) = instantiate_parameters(ch, reversals, conductances)
    eqs = [D(mH) ~ (1 / τm(ch, V)) * (m∞(ch, V) - mH),
        I ~ g * mH * (E - V)]
    current = [eqs[2]]
    defaultmap = [mH => ch.mH, g => ch.gH]
    return ComponentSystem(ch, ODESystem(eqs, t, [V, mH, I], [g, E]; observed=current, defaults=defaultmap, name=get_name(ch)))
end

#################### Leak current #########################

struct leak{D<:Real} <: FlowChannel(Leak)
    gLeak::D
end

function (ch::leak)()
    @variables V(t)
    (I,) = instantiate_variables(ch, currents)
    (E,), (g,) = instantiate_parameters(ch, reversals, conductances)
    eqs = [I ~ g * (E - V)]
    current = [eqs[1]]
    defaultmap = [g => ch.gLeak]
    return ComponentSystem(ch, ODESystem(eqs, t, [V, I], [g, E]; observed=current, defaults=defaultmap, name=get_name(ch)))
end


