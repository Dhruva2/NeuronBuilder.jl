"""
Plasticity rule can keep the original channel. 
Should be a component so it can be called as such from build_neuron
Does the needful for each regulated.
Stores extra plasticity-related variables separately
Callin a PlasticityRule() can then do the same as RegIonChannel, but using pr.channel, i.e. the channel stuff is separate from the plasticity rule in the namespace
"""


"""
Extra sensors: calcium
extra parameters: τmRNa
replace conductance parameter with a state
give that state dynamics:

Differential(t)(gNa(t)) ~ (mRna_Na(t) - gNa(t)) / τg
Differential(t)(mRna_Na(t)) ~ (Ca_tgt - Ca(t)) / τNa

Auxiliary: makes mrna name: channel name + mRna
"""
struct OLearyCalcRegulation{T} <: PlasticityRule{Calcium}
    τmRNA::T
    τg::T
    Ca_tgt::T
end

function OLeary_reg(ch, τmRNA, τg, Ca_tgt) 
    !(isLeak(ch)) ? PlasticisedChannel(ch, OLearyCalcRegulation(τmRNA, τg, Ca_tgt)) : ch
end

"""
make new componentsystem from old component system.
- add sensors, e,g, Calcium
- make an augmented name like cas_regulated <= cas.
- flatten the namespace so that there is just one system, not plasticityrulesystem.channel_system
- change the dynamics trying to access and set variables abstractly
"""

function (o::OLearyCalcRegulation)(ch::FlowChannel)
    sys = ch().sys
    _eqs, _states, _parameters, _defaults, _observed = map(
        (equations, states, parameters, ModelingToolkit.get_defaults, observed)) do f
        f(sys)
    end

    old_conductance = getproperty(sys, conductances(ch)...; namespace=false)
    (var_conductance,) = instantiate_variables(ch, conductances)

    # potential bug to fix: only works for a single element of conducts
    _param_indxs = findall(x -> isequal(old_conductance, x), _parameters)

    ## turn conductances into variables
    substitution = Dict(_parameters[_param_indxs]... .=> var_conductance...)
    new_eqs = [substitute(eq, substitution) for eq in _eqs]

    plast_ch = PlasticisedChannel(ch, o)
    #order matters
    (sensed(plast_ch)[1] == Calcium) ? (Ca, E) = instantiate_variables(plast_ch, sensedvars) : (E, Ca) = instantiate_variables(plast_ch, sensedvars)
    @variables mRNA(t)
    τmRNA, Ca_tgt, τg = get_parameters(o)

    extended_eqs = [D(mRNA) ~ (Ca_tgt - Ca) / τmRNA,
        D(var_conductance) ~ (mRNA - var_conductance) / τg]
    new_eqs = vcat(
        new_eqs, extended_eqs
    )
    new_states = add_states(_states, [Ca, var_conductance, mRNA])

    new_params = vcat(deleteat!(_parameters,_param_indxs), [τmRNA, Ca_tgt, τg])

    _defaults[var_conductance] = pop!(_defaults, old_conductance)
    new_defaultmap = merge(_defaults, Dict(mRNA => 0.0), default_params(o))

    plast_sys = ODESystem(new_eqs, t, new_states, new_params; observed=_observed, defaults=new_defaultmap, name=Symbol(get_name(ch), :_regul))
    return ComponentSystem(plast_ch, plast_sys)
end


add_states(_states::Vector, add_me::Vector) = vcat(_states, add_me)

get_parameters(::OLearyCalcRegulation) = @parameters τmRNA, Ca_tgt, τg
default_params(o::OLearyCalcRegulation) = Dict(get_parameters(o) .=> (o.τmRNA, o.Ca_tgt, o.τg))

function (p::PlasticisedChannel)()
    p.mutation(p.channel)
end

#####################################################################
struct FranciCalcRegulation{T,M,V} <: PlasticityRule{Calcium}
    τmRNA::T
    τg::T
    Ca_tgt::T
    A::M
    channels::V
end

Franci_reg(ch, τmRNA, τg, Ca_tgt, A, v) = PlasticisedChannel(ch, FranciCalcRegulation(τmRNA, τg, Ca_tgt, A, v))

"""
Only need mRNAs from other channels, which are created in the plasticised channels
"""
function (f::FranciCalcRegulation)(ch::FlowChannel)

    sys = ch().sys
    _eqs, _states, _parameters, _defaults, _observed = map(
        (equations, states, parameters, ModelingToolkit.get_defaults, observed)) do f
        f(sys)
    end

    old_conductance = getproperty(sys, conductances(ch)...; namespace=false)
    (var_conductance,) = instantiate_variables(ch, conductances)

    v = f.channels
    i = ch_index(ch, v)
    vec_mRNAs = instantiate_variables(Symbol.(:mRNA_, get_name.(v)))
    mRNA = vec_mRNAs[i]

    # potential bug to fix: only works for a single element of conducts
    _param_indxs = findall(x -> isequal(old_conductance, x), _parameters)

    ## turn conductances into variables
    substitution = Dict(_parameters[_param_indxs]... .=> var_conductance...)
    new_eqs = [substitute(eq, substitution) for eq in _eqs]

    plast_ch = PlasticisedChannel(ch, f)
    (E, Ca) = instantiate_variables(plast_ch, sensedvars)

    τmRNA, Ca_tgt, τg, _ = get_parameters(f)
    # getting 8x8 = 64 parameters a_ij is too complicated so we take the values of A directly in the equations
    #only need A[i,:] row 

    extended_eqs = [D(mRNA) ~ (Ca_tgt - Ca + f.A[i,:]' * vec_mRNAs) / τmRNA, #d_y, d_m_i = rand, rand
        D(var_conductance) ~ (mRNA - var_conductance + rand()) / τg] #d_g_i = rand
    new_eqs = vcat(
        new_eqs, extended_eqs
    )
    new_states = add_states(_states, [var_conductance, vec_mRNAs...])

    new_params = [τmRNA, Ca_tgt, τg] #no A !

    _defaults[var_conductance] = pop!(_defaults, old_conductance)
    new_defaultmap = merge(_defaults, Dict(mRNA => 0.0), default_params(f))

    plast_sys = ODESystem(new_eqs, t, new_states, new_params; observed=_observed, defaults=new_defaultmap, name=Symbol(get_name(ch), :_networegul))
    return ComponentSystem(plast_ch, plast_sys)
end

function ch_index(ch::FlowChannel,v::Vector{FlowChannel})
    findfirst(x->x==ch,v)
end

get_parameters(::FranciCalcRegulation) = @parameters τmRNA, Ca_tgt, τg
default_params(f::FranciCalcRegulation) = Dict(get_parameters(f) .=> (f.τmRNA, f.Ca_tgt, f.τg))
