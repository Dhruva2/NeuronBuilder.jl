# struct regu_dyn
#     conductances::AbstractVector
#     A::Matrix
# end

abstract type Model end
struct Zerodominant <: Model end #trait that highlights the regulation mechanism of a type 
struct Onedominant <: Model end

#trait functions 
model(::Type{<:IonChannel}) = Zerodominant()
model(::Type{<:Array{IonChannel}}) = Onedominant()

#trait dispatch
regul_dyn(xs::T, Ca) where {T} = regul_dyn(model(T), xs, Ca)
regul_dyn(::Zerodominant, xs, Ca) = regul_indep(xs, Ca)
regul_dyn(::Onedominant, xs, Ca) = regul_network(xs, Ca)


function regul_indep(ch::IonChannel, Ca)
    g = fieldnames(typeof(ch))[1]
    R = fieldnames(typeof(ch))[2]
    τ_ch = external_params(ch)[2]
    parameters = @parameters τg Ca_tgt $τ_ch
    states = @variables $g(t) $R(t)
    conductance, mRna = states[1], states[2]
    τch = parameters[3]
    eqs = [D(mRna) ~ (Ca_tgt - Ca) / τch,
        D(conductance) ~ (mRna - conductance) / τg]
    #observed = g ?
    defaultmap = [conductance => getg(ch), mRna => getR(ch)]
    return eqs, states, parameters, defaultmap
end

#function regul_network(conductances::RegNetwork, Ca)