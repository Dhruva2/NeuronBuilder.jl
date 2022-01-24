# struct regu_dyn
#     conductances::AbstractVector
#     A::Matrix
# end

abstract type Model end
struct Zerodominant <: Model end #trait that highlights the regulation mechanism of a type 
struct Onedominant <: Model end

#trait functions 
model(::Type{<:RegIonChannel}) = Zerodominant()
model(::Type{<:Array{IonChannel}}) = Onedominant()

#trait dispatch
regul_dyn(xs::T, Ca) where {T} = regul_dyn(model(T), xs, Ca)
regul_dyn(::Zerodominant, xs, Ca) = regul_indep(xs, Ca)
regul_dyn(::Onedominant, xs, Ca) = regul_network(xs, Ca)

function regul_indep(ch::RegIonChannel, Ca)
    g = get_g(ch.ionch)
    R = Symbol(:R,get_name(ch.ionch))
    τ_ch = ch.τion
    parameters = @parameters $τ_ch Ca_tgt τg
    states = @variables $g(t) $R(t)
    conductance, mRna = states[1], states[2]
    τch = parameters[1]
    eqs = [D(mRna) ~ (Ca_tgt - Ca) / τch,
        D(conductance) ~ (mRna - conductance) / τg]
    defaultmap = [conductance => getfield(ch.ionch,get_g(ch.ionch)), mRna => ch.Rion]
    return eqs, states, parameters, defaultmap
end

#function regul_network(conductances::RegNetwork, Ca)