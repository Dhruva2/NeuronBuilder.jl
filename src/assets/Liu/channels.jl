mutable struct RegIon <: RegIonChannel
    ionch::IonChannel
end

getR(ch) = getfield(ch, fieldnames(typeof(ch))[2])

#################### NaV ###############################
mutable struct NaV{S,T} <: IonChannel
    gNa::S # scalar or variable
    RNa::S
    mNa::T
    hNa::T
end

NaV(x; reg = false) = reg ? RegIon(NaV(x, x, 0.0, 0.0)) : NaV(x, 0.0, 0.0, 0.0) #if regulated, R is initialized with the same value as g
ionic_current(::NaV, sys::ODESystem) = sys.INa
external_params(::NaV) = (:ENa, :τNa, :Ca_tgt, :τg)
m∞(::NaV, V) = 1.0 / (1.0 + exp((V + 25.5) / -5.29))
h∞(::NaV, V) = 1.0 / (1.0 + exp((V + 48.9) / 5.18))
τm(::NaV, V) = 1.32 - 1.26 / (1 + exp((V + 120.0) / -25.0))
τh(::NaV, V) = (0.67 / (1.0 + exp((V + 62.9) / -10.0))) * (1.5 + 1.0 / (1.0 + exp((V + 34.9) / 3.6)))

function channel_dynamics(ch::NaV, V, Ca)
    if getR(ch) == 0.0
        states = @variables mNa(t) hNa(t) INa(t)
        parameters = @parameters gNa ENa
    elseif getR(ch) != 0.0
        states = @variables mNa(t) hNa(t) INa(t) gNa(t)
        parameters = @parameters ENa
    end
    eqs = [D(mNa) ~ (1 / τm(ch, V)) * (m∞(ch, V) - mNa),
        D(hNa) ~ (1 / τh(ch, V)) * (h∞(ch, V) - hNa),
        INa ~ gNa * mNa^3 * hNa * (ENa - V)]
    current = [eqs[3]]
    defaultmap = [mNa => ch.mNa, hNa => ch.hNa, gNa => ch.gNa]
    return eqs, states, parameters, current, defaultmap
end

#################### Slow calcium current #############################

mutable struct CaS{S,T} <: IonChannel
    gCaS::S
    RCaS::S
    mCaS::T
    hCaS::T
end

CaS(x; reg = false) = reg ? RegIon(CaS(x, x, 20.0, 0.0)) : CaS(x, 0.0, 20.0, 0.0)
CaS(x, y; reg = false) = reg ? RegIon(CaS(x, x, y, 0.0)) : CaS(x, 0.0, y, 0.0)
ECa(::CaS, Ca) = (500.0) * (8.6174e-5) * (283.15) * (log(max((3000.0 / Ca), 0.001)))
external_params(::CaS) = (nothing, :τCaS, :Ca_tgt, :τg)
ionic_current(::CaS, sys::ODESystem) = sys.ICaS
calcium_current(::CaS, sys::ODESystem) = sys.ICaS
m∞(::CaS, V) = 1.0 / (1.0 + exp((V + 33.0) / -8.1))
h∞(::CaS, V) = 1.0 / (1.0 + exp((V + 60.0) / 6.2))
τm(::CaS, V) = 1.4 + 7.0 / (exp((V + 27.0) / 10.0) + exp((V + 70.0) / -13.0))
τh(::CaS, V) = 60.0 + 150.0 / (exp((V + 55.0) / 9.0) + exp((V + 65.0) / -16.0))

function channel_dynamics(ch::CaS, V, Ca)
    if getR(ch) == 0.0
        states = @variables mCaS(t) hCaS(t) ICaS(t)
        parameters = @parameters gCaS
    elseif getR(ch) == 0.0
        states = @variables mCaS(t) hCaS(t) ICaS(t) gCaS(t)
        parameters = ()
    end
    eqs = [D(mCaS) ~ (1 / τm(ch, V)) * (m∞(ch, V) - mCaS),
        D(hCaS) ~ (1 / τh(ch, V)) * (h∞(ch, V) - hCaS),
        ICaS ~ gCaS * mCaS^3 * hCaS * (ECa(ch, Ca) - V)]
    current = [eqs[3]]
    defaultmap = [mCaS => ch.mCaS, hCaS => ch.hCaS, gCaS => ch.gCaS]
    return eqs, states, parameters, current, defaultmap
end

#################### Transient calcium current ######################

mutable struct CaT{S,T} <: IonChannel
    gCaT::S
    RCaT::S
    mCaT::T
    hCaT::T
end

CaT(x; reg = false) = reg ? RegIon(CaT(x, x, 0.0, 0.0)) : CaT(x, 0.0, 0.0, 0.0)
ionic_current(::CaT, sys::ODESystem) = sys.ICaT
calcium_current(::CaT, sys::ODESystem) = sys.ICaT
ECa(::CaT, Ca) = (500.0) * (8.6174e-5) * (283.15) * (log(max((3000.0 / Ca), 0.001)))
external_params(::CaT) = (nothing, :τCaT, :Ca_tgt, :τg)
m∞(::CaT, V) = 1.0 / (1.0 + exp((V + 27.1) / -7.2))
h∞(::CaT, V) = 1.0 / (1.0 + exp((V + 32.1) / 5.5))
τm(::CaT, V) = 21.7 - 21.3 / (1.0 + exp((V + 68.1) / -20.5));
τh(::CaT, V) = 105.0 - 89.8 / (1.0 + exp((V + 55.0) / -16.9));

function channel_dynamics(ch::CaT, V, Ca)
    if getR(ch) == 0.0
        states = @variables mCaT(t) hCaT(t) ICaT(t)
        parameters = @parameters gCaT
    elseif getR(ch) != 0.0
        states = @variables mCaT(t) hCaT(t) ICaT(t) gCaT(t)
        parameters = @parameters ()
    end
    eqs = [D(mCaT) ~ (1 / τm(ch, V)) * (m∞(ch, V) - mCaT),
        D(hCaT) ~ (1 / τh(ch, V)) * (h∞(ch, V) - hCaT),
        ICaT ~ gCaT * mCaT^3 * hCaT * (ECa(ch, Ca) - V)]
    current = [eqs[3]]
    defaultmap = [mCaT => ch.mCaT, hCaT => ch.hCaT, gCaT => ch.gCaT]
    return eqs, states, parameters, current, defaultmap
end

#################### A-type potassium current #########################

mutable struct Ka{S,T} <: IonChannel
    gKa::S #can be a scalar or a variable
    RKa::S #if unregulated then this is 0.0
    mKa::T
    hKa::T
end

Ka(x; reg = false) = reg ? RegIon(Ka(x, x, 0.0, 0.0)) : Ka(x, 0.0, 0.0, 0.0)
ionic_current(::Ka, sys::ODESystem) = sys.IKa
external_params(::Ka) = (:EK, :τKa, :Ca_tgt, :τg)
m∞(::Ka, V) = 1.0 / (1.0 + exp((V + 27.2) / -8.7))
h∞(::Ka, V) = 1.0 / (1.0 + exp((V + 56.9) / 4.9))
τm(::Ka, V) = 11.6 - 10.4 / (1.0 + exp((V + 32.9) / -15.2))
τh(::Ka, V) = 38.6 - 29.2 / (1.0 + exp((V + 38.9) / -26.5))

function channel_dynamics(ch::Ka, V, Ca)
    if getR(ch) == 0.0
        states = @variables mKa(t) hKa(t) IKa(t)
        parameters = @parameters gKa EK
    elseif getR(ch) != 0.0
        states = @variables mKa(t) hKa(t) IKa(t) gKa(t)
        parameters = @parameters EK
    end
    eqs = [D(mKa) ~ (1 / τm(ch, V)) * (m∞(ch, V) - mKa),
        D(hKa) ~ (1 / τh(ch, V)) * (h∞(ch, V) - hKa),
        IKa ~ gKa * mKa^3 * hKa * (EK - V)]
    current = [eqs[3]]
    defaultmap = [mKa => ch.mKa, hKa => ch.hKa, gKa => ch.gKa]
    return eqs, states, parameters, current, defaultmap
end

################### Calcium-activated potassium current ########

mutable struct KCa{S,T} <: IonChannel
    gKCa::S
    RKCa::S
    mKCa::T
end

KCa(x; reg = false) = reg ? RegIon(KCa(x, x, 0.0)) : KCa(x, 0.0, 0.0)
ionic_current(::KCa, sys::ODESystem) = sys.IKCa
external_params(::KCa) = (:EK, :τKCa, :Ca_tgt, :τg)
m∞(::KCa, V, Ca) = (Ca / (Ca + 3.0)) / (1.0 + exp((V + 28.3) / -12.6));
τm(::KCa, V) = 90.3 - 75.1 / (1.0 + exp((V + 46.0) / -22.7));

function channel_dynamics(ch::KCa, V, Ca)
    if getR(ch) == 0.0
        states = @variables mKCa(t) IKCa(t)
        parameters = @parameters gKCa EK
    elseif getR(ch) != 0.0
        states = @variables mKCa(t) IKCa(t) gKCa(t)
        parameters = @parameters EK
    end
    eqs = [D(mKCa) ~ (1 / τm(ch, V)) * (m∞(ch, V, Ca) - mKCa),
        IKCa ~ (gKCa * mKCa^4) * (EK - V)]
    current = [eqs[2]]
    defaultmap = [mKCa => ch.mKCa, gKCa => ch.gKCa]
    return eqs, states, parameters, current, defaultmap
end

#################### Delayed rectifier potassium current ######################
mutable struct Kdr{S,T} <: IonChannel
    gKdr::S
    RKdr::S
    mKdr::T
end

Kdr(x; reg = false) = reg ? RegIon(Kdr(x, x, 0.0)) : Kdr(x, 0.0, 0.0)
ionic_current(::Kdr, sys::ODESystem) = sys.IKdr
external_params(::Kdr) = (:EK, :τKdr, :Ca_tgt, :τg)
m∞(::Kdr, V) = 1.0 / (1.0 + exp((V + 12.3) / -11.8));
τm(::Kdr, V) = 7.2 - 6.4 / (1.0 + exp((V + 28.3) / -19.2));

function channel_dynamics(ch::Kdr, V, Ca)
    if getR(ch) == 0.0
        states = @variables mKdr(t) IKdr(t)
        parameters = @parameters gKdr EK
    elseif getR(ch) != 0.0
        states = @variables mKdr(t) IKdr(t) gKdr(t)
        parameters = @parameters EK
    end
    eqs = [D(mKdr) ~ (1 / τm(ch, V)) * (m∞(ch, V) - mKdr),
        IKdr ~ (gKdr * mKdr^4) * (EK - V)]
    current = [eqs[2]]
    defaultmap = [mKdr => ch.mKdr, gKdr => ch.gKdr]
    return eqs, states, parameters, current, defaultmap
end

#################### H current ####################

mutable struct H{S,T} <: IonChannel
    gH::S
    RH::S
    mH::T
end

H(x; reg = false) = reg ? RegIon(H(x, x, 0.0)) : H(x, 0.0, 0.0)
ionic_current(::H, sys::ODESystem) = sys.IH
external_params(::H) = (:EH, :τH, :Ca_tgt, :τg)
m∞(::H, V) = 1.0 / (1.0 + exp((V + 70.0) / 6.0))
τm(::H, V) = (272.0 + 1499.0 / (1.0 + exp((V + 42.2) / -8.73)))

function channel_dynamics(ch::H, V, Ca)
    if getR(ch) == 0.0
        states = @variables mH(t) IH(t)
        parameters = @parameters gH EH
    elseif getR(ch) != 0.0
        states = @variables mH(t) IH(t) gH(t)
        parameters = @parameters EH
    end
    eqs = [D(mH) ~ (1 / τm(ch, V)) * (m∞(ch, V) - mH),
        IH ~ gH * mH * (EH - V)]
    current = [eqs[2]]
    defaultmap = [mH => ch.mH, gH => ch.gH]
    return eqs, states, parameters, current, defaultmap
end


#################### Leak current #########################

mutable struct Leak{Float64} <: IonChannel
    gLeak::Float64
end

ionic_current(::Leak, sys::ODESystem) = sys.ILeak
external_params(::Leak) = (:ELeak,)

function channel_dynamics(ch::Leak, V, Ca)
    states = @variables ILeak(t)
    parameters = @parameters gLeak ELeak
    eqs = [ILeak ~ gLeak * (ELeak - V)]
    current = [eqs[1]]
    defaultmap = [gLeak => ch.gLeak]
    return eqs, states, parameters, current, defaultmap
end
