abstract type Model end
struct Zerodominant <: Model end #trait that highlights the regulation mechanism of a type 
# struct Onedominant <: Model end
## future work: regulation network for conductances. Onedominant is not used for now.

## trait functions 
model(::Type{<:RegIonChannel}) = Zerodominant()
#model(::Type{<:Array{IonChannel}}) = Onedominant()

##trait dispatch
regul_dyn(xs::T, Ca) where {T} = regul_dyn(model(T), xs, Ca)
regul_dyn(::Zerodominant, xs, Ca) = regul_indep(xs, Ca)
#regul_dyn(::Onedominant, xs, Ca) = regul_network(xs, Ca)

function regul_indep(ch::RegIonChannel, Ca)
    g = get_g(ch.ionch)
    R = get_mRna_name(ch)
    τ_ch = ch.τion
    parameters = @parameters $τ_ch Ca_tgt τg
    states = @variables $g(t) $R(t)
    conductance, mRna = states[1], states[2]
    τch = parameters[1]
    eqs = [D(mRna) ~ (Ca_tgt - Ca) / τch,
        D(conductance) ~ (mRna - conductance) / τg]
    defaultmap = [conductance => getfield(ch.ionch, get_g(ch.ionch)), mRna => ch.Rion]
    return eqs, states, parameters, defaultmap
end

#function regul_network(conductances::RegNetwork, Ca)

function Regulated(i::IonChannel)
    if get_name(i) == :Leak
        return i
    else
        R = getfield(i, get_g(i)) #set initial value for R same as g
        return NeuronBuilder.RegIon(i, R, Symbol(:τ, get_name(i)))
    end
end

external_params(ch::RegIonChannel) = (external_params(ch.ionch)..., :Ca_tgt, :τg)
get_mRna_name(ch::RegIonChannel) = Symbol(:mRna, :_, get_name(ch.ionch))
ionic_current(ch::RegIonChannel, sys::ODESystem) = ionic_current(ch.ionch, sys)
calcium_current(ch::RegIonChannel, sys::ODESystem) = calcium_current(ch.ionch, sys)