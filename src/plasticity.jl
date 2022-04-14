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

OLeary_reg(ch, τmRNA, τg, Ca_tgt) = PlasticisedChannel(ch, OLearyCalcRegulation(τmRNA, τg, Ca_tgt))

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

    old_conductance = getproperty(sys, conductance(ch)...; namespace=false)
    (var_conductance,) = instantiate_variables(ch, conductance)

    # potential bug to fix: only works for a single element of conducts
    _param_indxs = findall(x -> isequal(old_conductance, x), _parameters)

    ## turn conductances into variables
    substitution = Dict(_parameters[_param_indxs]... .=> var_conductance...)
    new_eqs = [substitute(eq, substitution) for eq in _eqs]

    plast_ch = PlasticisedChannel(ch, o)
    (E, Ca) = instantiate_variables(plast_ch, sensedvars)
    @variables mRNA(t)
    τmRNA, Ca_tgt, τg = get_parameters(o)

    extended_eqs = [D(mRNA) ~ (Ca_tgt - Ca) / τmRNA,
        D(var_conductance) ~ (mRNA - var_conductance) / τg]
    new_eqs = vcat(
        new_eqs, extended_eqs
    )
    new_states = add_states(_states, var_conductance)

    new_params = [τmRNA, Ca_tgt, τg]

    _defaults[var_conductance] = pop!(_defaults, old_conductance)
    new_defaultmap = merge(_defaults, Dict(mRNA => 0.0), default_params(o))

    plast_sys = ODESystem(new_eqs, t, new_states, new_params; observed=_observed, defaults=new_defaultmap, name=Symbol(get_name(ch), :_regul))
    return ComponentSystem(plast_ch, plast_sys)
end


add_states(_states, add_me) = vcat(_states, add_me)

get_parameters(::OLearyCalcRegulation) = @parameters τmRNA, Ca_tgt, τg
default_params(o::OLearyCalcRegulation) = Dict(get_parameters(o) .=> (o.τmRNA, o.Ca_tgt, o.τg))


function (p::PlasticisedChannel)()
    p.mutation(p.channel)
end