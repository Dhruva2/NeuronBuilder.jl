struct Regulated end 

struct regu_dyn
    conductances::AbstractVector
    A::Matrix
end

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
    τch = external_params(ch)[2]
    parameters = @parameters τg Ca_tgt $τch
    states = @variables $g(t) $R(t)
    conductance,mRna = states[1], states[2]
    eqs = [D(mRna) ~ (Ca_tgt - Ca) / τch,
        D(conductance) ~ (mRna - conductance) / τg]
    #observed = g ?
    defaultmap = [conductance => getg(ch), mRna => getR(ch)]
    @named gandR = ODESystem(eqs, t, states, parameters, defaults = defaultmap)
    return gandR
end

#function regul_network(conductances::RegNetwork, Ca)

ionic_current(ion_chan::Channel) = ion_chan.sys.states[end]

calcium_current(IonChannel) = hasCa(ion) ? ionic_current(ch) : Num(0)

#calcium_current()