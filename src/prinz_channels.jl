#################### NaV ###############################
mutable struct Prinz_NaV{S,T} <: IonChannel
    ḡNa::S
    mNa::T
    hNa::T
    #ENa::T
end

Prinz_NaV(x) = Prinz_NaV(x, 0.,0.)


m∞(::Prinz_NaV, V) =  1.0/(1.0+exp((V+25.5)/-5.29))
h∞(::Prinz_NaV, V) =  1.0/(1.0+exp((V+48.9)/5.18))
τm(::Prinz_NaV, V) =  1.32 - 1.26/(1+exp((V+120.0)/-25.0))
τh(::Prinz_NaV, V) =  (1.34/(1.0+exp((V+62.9)/-10.0)))*(1.5+1.0/(1.0+exp((V+34.9)/3.6)))
ionic_current(::Prinz_NaV, sys::ODESystem) = sys.INa
external_params(::Prinz_NaV) = (:ENa,)

function channel_dynamics(ch::Prinz_NaV, V, Ca, D, t)    
    states = @variables mNa(t) hNa(t) INa(t)
    parameters = @parameters ḡNa ENa 
    eqs = [ D(mNa) ~       (1/τm(ch, V))*(m∞(ch, V) - mNa), 
            D(hNa) ~       (1/τh(ch, V))*(h∞(ch, V) - hNa),
            INa ~ ḡNa*mNa^3*hNa*(ENa - V)]   
    current = [eqs[3]]
    defaultmap = [mNa => ch.mNa, hNa => ch.hNa, ḡNa => ch.ḡNa]
    return eqs, states, parameters, current, defaultmap
end  


#################### Slow calcium current #############################
mutable struct Prinz_CaS{S,T} <: IonChannel
    ḡCaS::S
    mCaS::T
    hCaS::T
end

Prinz_CaS(x,y) = Prinz_CaS(x,y, 0.)
Prinz_CaS(x) = Prinz_CaS(x,20., 0.)


m∞(::Prinz_CaS, V) = 1.0/(1.0+exp((V+33.0)/-8.1))
h∞(::Prinz_CaS, V) = 1.0/(1.0+exp((V+60.0)/6.2))
τm(::Prinz_CaS, V) = 2.8 + 14.0/(exp((V+27.0)/10.0) + exp((V+70.0)/-13.0))
τh(::Prinz_CaS, V) = 120.0 + 300.0/(exp((V+55.0)/9.0) + exp((V+65.0)/-16.0))

ionic_current(::Prinz_CaS, sys::ODESystem) = sys.ICaS
calcium_current(::Prinz_CaS, sys::ODESystem) = sys.ICaS
ECa(::Prinz_CaS, Ca) = (500.0)*(8.6174e-5)*(283.15)*(log(max((3000.0/Ca), 0.001)))
external_params(::Prinz_CaS) = nothing

function channel_dynamics(ch::Prinz_CaS, V, Ca, D, t)
    states = @variables mCaS(t) hCaS(t) ICaS(t)
    parameters = @parameters ḡCaS
    eqs = [ D(mCaS) ~       (1/τm(ch, V))*(m∞(ch, V) - mCaS), 
            D(hCaS) ~       (1/τh(ch, V))*(h∞(ch, V) - hCaS),
            ICaS    ~       ḡCaS*mCaS^3*hCaS*(ECa(ch, Ca) - V)
            ]
    current = [eqs[3]]
    defaultmap = [mCaS => ch.mCaS, hCaS => ch.hCaS, ḡCaS => ch.ḡCaS]
    return eqs, states, parameters, current, defaultmap
end

#################### Transient calcium current ######################

mutable struct Prinz_CaT{S,T} <: IonChannel
    ḡCaT::S
    mCaT::T
    hCaT::T
end
Prinz_CaT(x) = Prinz_CaT(x, 0., 0.)
m∞(::Prinz_CaT, V) = 1.0/(1.0 + exp((V+27.1)/-7.2))
h∞(::Prinz_CaT, V) = 1.0/(1.0 + exp((V+32.1)/5.5))
τm(::Prinz_CaT, V) = 43.4 - 42.6/(1.0 + exp((V+68.1)/-20.5))
τh(::Prinz_CaT, V) = 210.0 - 179.6/(1.0 + exp((V+55.0)/-16.9))


ionic_current(::Prinz_CaT, sys::ODESystem) = sys.ICaT
calcium_current(::Prinz_CaT, sys::ODESystem) = sys.ICaT
ECa(::Prinz_CaT, Ca) = (500.0)*(8.6174e-5)*(283.15)*(log(max((3000.0/Ca), 0.001)))
external_params(::Prinz_CaT) = nothing

function channel_dynamics(ch::Prinz_CaT, V, Ca, D, t)
    states = @variables mCaT(t) hCaT(t) ICaT(t)
    parameters = @parameters ḡCaT
    eqs = [ D(mCaT) ~       (1/τm(ch, V))*(m∞(ch, V) - mCaT), 
            D(hCaT) ~       (1/τh(ch, V))*(h∞(ch, V) - hCaT),
            ICaT ~ ḡCaT*mCaT^3*hCaT*(ECa(ch, Ca) - V),
            ]
    current = [eqs[3]]
    defaultmap = [mCaT => ch.mCaT, hCaT => ch.hCaT, ḡCaT => ch.ḡCaT]
    return eqs, states, parameters, current, defaultmap
end

####################  #########################
"""
A-type potassium current
"""
mutable struct Prinz_Ka{S,T} <: IonChannel
    ḡKa::S
    mKa::T
    hKa::T
    #EK::T -80
end
Prinz_Ka(x) = Prinz_Ka(x,0.,0.)


m∞(::Prinz_Ka, V) = 1.0/(1.0+exp((V+27.2)/-8.7))
h∞(::Prinz_Ka, V) = 1.0/(1.0+exp((V+56.9)/4.9))
τm(::Prinz_Ka, V) = 23.2 - 20.8/(1.0+exp((V+32.9)/-15.2))
τh(::Prinz_Ka, V) = 77.2 - 58.4/(1.0+exp((V+38.9)/-26.5))

ionic_current(::Prinz_Ka, sys::ODESystem) = sys.IKa
external_params(::Prinz_Ka) = (:EK,)

function channel_dynamics(ch::Prinz_Ka, V, Ca, D, t)    
    states = @variables mKa(t) hKa(t) IKa(t)
    parameters = @parameters ḡKa EK 
    eqs = [ D(mKa) ~       (1/τm(ch, V))*(m∞(ch, V) - mKa), 
            D(hKa) ~       (1/τh(ch, V))*(h∞(ch, V) - hKa),
            IKa ~ ḡKa*mKa^3*hKa*(EK - V)]   
    current = [eqs[3]]
    defaultmap = [mKa => ch.mKa, hKa => ch.hKa, ḡKa => ch.ḡKa]
    return eqs, states, parameters, current, defaultmap
end  

################### Calcium-activated potassium current ########

mutable struct Prinz_KCa{S,T} <: IonChannel
    ḡKCa::S
    mKCa::T
    #EK::T
end

Prinz_KCa(x) = Prinz_KCa(x,0.)
m∞(::Prinz_KCa, V, Ca) = (Ca/(Ca+3.0))/(1.0+exp((V+28.3)/-12.6));
τm(::Prinz_KCa, V) = 180.6 - 150.2/(1.0+exp((V+46.0)/-22.7))
ionic_current(::Prinz_KCa, sys::ODESystem) = sys.IKCa
external_params(::Prinz_KCa) = (:EK,)

function channel_dynamics(ch::Prinz_KCa, V, Ca, D, t)    
    states = @variables mKCa(t) IKCa(t)
    parameters = @parameters ḡKCa EK 
    eqs = [ D(mKCa) ~       (1/τm(ch, V))*(m∞(ch, V, Ca) - mKCa), 
            IKCa ~ (ḡKCa*mKCa^4)*(EK - V)]   
    current = [eqs[2]]
    defaultmap = [mKCa => ch.mKCa, ḡKCa => ch.ḡKCa]
    return eqs, states, parameters, current, defaultmap
end  



"""
    Delayed rectifier potassium current
"""
mutable struct Prinz_Kdr{S,T} <: IonChannel
    ḡKdr::S
    mKdr::T
    #EK::T
end
Prinz_Kdr(x) = Prinz_Kdr(x,0.)


m∞(::Prinz_Kdr, V)=  1.0/(1.0+exp((V+12.3)/-11.8));
τm(::Prinz_Kdr, V)=  14.4 - 12.8/(1.0+exp((V+28.3)/-19.2))
ionic_current(::Prinz_Kdr, sys::ODESystem) = sys.IKdr
external_params(::Prinz_Kdr) = (:EK,)


function channel_dynamics(ch::Prinz_Kdr, V, Ca, D, t)    
    states = @variables mKdr(t) IKdr(t)
    parameters = @parameters ḡKdr EK 
    eqs = [ D(mKdr) ~       (1/τm(ch, V))*(m∞(ch, V) - mKdr), 
            IKdr ~ (ḡKdr*mKdr^4)*(EK - V)]   
    current = [eqs[2]]
    defaultmap = [mKdr => ch.mKdr, ḡKdr => ch.ḡKdr]
    return eqs, states, parameters, current, defaultmap
end  

"""
H current
"""
mutable struct Prinz_H{S,T} <: IonChannel
    ḡH::S
    mH::T
    # EH::T -20.
end

Prinz_H(x) = Prinz_H(x, 0.)


m∞(::Prinz_H, V) = 1.0/(1.0+exp((V+75.0)/5.5))
τm(::Prinz_H, V) = (2/( exp((V+169.7)/(-11.6)) + exp((V- 26.7)/(14.3)) ))
ionic_current(::Prinz_H, sys::ODESystem) = sys.IH
external_params(::Prinz_H) = (:EH,)

function channel_dynamics(ch::Prinz_H, V, Ca, D, t)    
    states = @variables mH(t) IH(t)
    parameters = @parameters ḡH EH 
    eqs = [ D(mH) ~       (1/τm(ch, V))*(m∞(ch, V) - mH), 
            IH ~ ḡH*mH*(EH - V)]   
    current = [eqs[2]]
    defaultmap = [mH => ch.mH, ḡH => ch.ḡH]
    return eqs, states, parameters, current, defaultmap
end  


"""
leak current
"""
mutable struct Prinz_Leak{S} <: IonChannel
    ḡLeak::S 
    #ELeak::S -50
end

ionic_current(::Prinz_Leak, sys::ODESystem) = sys.ILeak
external_params(::Prinz_Leak) = (:ELeak,)

function channel_dynamics(ch::Prinz_Leak, V, Ca, D, t)    
    states = @variables ILeak(t)
    parameters = @parameters ḡLeak ELeak 
    eqs = [ILeak ~ ḡLeak*(ELeak - V)]   
    current = [eqs[1]]
    defaultmap = [ḡLeak => ch.ḡLeak]
    return eqs, states, parameters, current, defaultmap
end  