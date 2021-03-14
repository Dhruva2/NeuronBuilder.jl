abstract type Channel end

"""
return unparameterised type as symbol
"""
function get_name(ch::T) where T<:Channel 
    Base.typename(ch |> typeof) |> Symbol
end

"""
fallback option for channels without a calcium current
"""
function calcium_current(ch::T, sys::ODESystem) where T<:Channel
    return Num(0.)
end

#################### Channels ###############################


#################### NaV ###############################
mutable struct NaV{S,T} <: Channel
    ḡNa::S
    mNa::T
    hNa::T
    ENa::T
end

NaV() = NaV(100., 0.,0., 50.)
NaV(x,y) = NaV(x, 0.,0., y)

m∞(::NaV, V) =  1.0/(1.0+exp((V+25.5)/-5.29))
h∞(::NaV, V) =  1.0/(1.0+exp((V+48.9)/5.18))
τm(::NaV, V) =  1.32 - 1.26/(1+exp((V+120.0)/-25.0))
τh(::NaV, V) =  (0.67/(1.0+exp((V+62.9)/-10.0)))*(1.5+1.0/(1.0+exp((V+34.9)/3.6)))
voltage_current(::NaV, sys::ODESystem) = sys.INa


function channel_dynamics(ch::NaV, V, Ca, D, t)    
    states = @variables mNa(t) hNa(t) INa(t)
    parameters = @parameters ḡNa ENa 
    eqs = [ D(mNa) ~       (1/τm(ch, V))*(m∞(ch, V) - mNa), 
            D(hNa) ~       (1/τh(ch, V))*(h∞(ch, V) - hNa),
            INa ~ ḡNa*mNa^3*hNa*(ENa - V)]   
    current = [eqs[3]]
    u0map = [mNa => ch.mNa, hNa => ch.hNa]
    pmap = [ḡNa => ch.ḡNa, ENa => ch.ENa]
    return eqs, states, parameters, current, u0map, pmap
end  


#################### Slow calcium current #############################
mutable struct CaS{S,T} <: Channel
    ḡCaS::S
    τCa::S
    mCaS::T
    hCaS::T
end

CaS(x,y) = CaS(x,y, 0.,0.)
CaS(x) = CaS(x,20., 0., 0.)

m∞(::CaS, V) = 1.0/(1.0+exp((V+33.0)/-8.1))
h∞(::CaS, V) = 1.0/(1.0+exp((V+60.0)/6.2))
τm(::CaS, V) = 1.4 + 7.0/(exp((V+27.0)/10.0) + exp((V+70.0)/-13.0));
τh(::CaS, V) = 60.0 + 150.0/(exp((V+55.0)/9.0) + exp((V+65.0)/-16.0));
voltage_current(::CaS, sys::ODESystem) = sys.ICaS
calcium_current(::CaS, sys::ODESystem) = sys.ICaS_Ca
ECa(::CaS, Ca) = (500.0)*(8.6174e-5)*(283.15)*(log((3000.0/Ca)))


function channel_dynamics(ch::CaS, V, Ca, D, t)
    states = @variables mCaS(t) hCaS(t) ICaS(t) ICaS_Ca(t)
    parameters = @parameters ḡCaS τCa
    eqs = [ D(mCaS) ~       (1/τm(ch, V))*(m∞(ch, V) - mCaS), 
            D(hCaS) ~       (1/τh(ch, V))*(h∞(ch, V) - hCaS),
            ICaS ~ ḡCaS*mCaS^3*hCaS*(ECa(ch, Ca) - V),
            ICaS_Ca ~ (1. / τCa)* (0.025 + 0.94ICaS - 0.5Ca)
            ]
    current = eqs[3:4]
    u0map = [mCaS => ch.mCaS, hCaS => ch.hCaS]
    pmap = [ḡCaS => ch.ḡCaS, τCa => ch.τCa]
    return eqs, states, parameters, current, u0map, pmap
end

#################### Transient calcium current ######################

mutable struct CaT{S,T} <: Channel
    ḡCaT::S
    τCa::S
    mCaT::T
    hCaT::T
end
CaT(x,y) = CaT(x,y,0.,0.)
CaT(x) = CaT(x, 20., 0., 0.)
m∞(::CaT, V) = 1.0/(1.0 + exp((V+27.1)/-7.2))
h∞(::CaT, V) = 1.0/(1.0 + exp((V+32.1)/5.5))
τm(::CaT, V) = 21.7 - 21.3/(1.0 + exp((V+68.1)/-20.5));
τh(::CaT, V) = 105.0 - 89.8/(1.0 + exp((V+55.0)/-16.9));

voltage_current(::CaT, sys::ODESystem) = sys.ICaT
calcium_current(::CaT, sys::ODESystem) = sys.ICaT_Ca
ECa(::CaT, Ca) = (500.0)*(8.6174e-5)*(283.15)*(log((3000.0/Ca)))

function channel_dynamics(ch::CaT, V, Ca, D, t)
    states = @variables mCaT(t) hCaT(t) ICaT(t) ICaT_Ca(t)
    parameters = @parameters ḡCaT τCa
    eqs = [ D(mCaT) ~       (1/τm(ch, V))*(m∞(ch, V) - mCaT), 
            D(hCaT) ~       (1/τh(ch, V))*(h∞(ch, V) - hCaT),
            ICaT ~ ḡCaT*mCaT^3*hCaT*(ECa(ch, Ca) - V),
            ICaT_Ca ~ (1. / τCa)* (0.025 + 0.94ICaT - 0.5Ca)
            ]
    current = eqs[3:4]
    u0map = [mCaT => ch.mCaT, hCaT => ch.hCaT]
    pmap = [ḡCaT => ch.ḡCaT, τCa => ch.τCa]
    return eqs, states, parameters, current, u0map, pmap
end

#################### A-type potassium current #########################

mutable struct Ka{S,T} <: Channel
    ḡKa::S
    mKa::T
    hKa::T
    EK::T
end
Ka(x,y) = Ka(x,0., 0.,y)
Ka(x) = Ka(x,0.,0.,80.)

m∞(::Ka, V) = 1.0/(1.0+exp((V+27.2)/-8.7))
h∞(::Ka, V) = 1.0/(1.0+exp((V+56.9)/4.9))
τm(::Ka, V) = 11.6 - 10.4/(1.0+exp((V+32.9)/-15.2));
τh(::Ka, V) = 38.6 - 29.2/(1.0+exp((V+38.9)/-26.5));
voltage_current(::Ka, sys::ODESystem) = sys.IKa

function channel_dynamics(ch::Ka, V, Ca, D, t)    
    states = @variables mKa(t) hKa(t) IKa(t)
    parameters = @parameters ḡKa EK 
    eqs = [ D(mKa) ~       (1/τm(ch, V))*(m∞(ch, V) - mKa), 
            D(hKa) ~       (1/τh(ch, V))*(h∞(ch, V) - hKa),
            IKa ~ ḡKa*mKa^3*hKa*(EK - V)]   
    current = [eqs[3]]
    u0map = [mKa => ch.mKa, hKa => ch.hKa]
    pmap = [ḡKa => ch.ḡKa, EK => ch.EK]
    return eqs, states, parameters, current, u0map, pmap
end  

################### Calcium-activated potassium current ########



#  cn = CalciumNeuron(50.,-20.,-80.,-50.,20.)
                        #ENa Eh EK Eleak τCa)