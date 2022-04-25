
function (p::PlasticisedChannel)(n::Neuron)
    p.modification(p.channel, n)
end

struct OLearyCalcRegulation{T} <: PlasticityRule{Calcium}
    τmRNA::T
    τg::T
    Ca_tgt::T
end

OLeary_reg(ch, τmRNA, τg, Ca_tgt) = PlasticisedChannel(ch, OLearyCalcRegulation(τmRNA, τg, Ca_tgt))::FlowChannel

"""
make new component system from old component system.
- add variables for new sensors, eg, Calcium
- param_to_state substitutes new variables for old parameters in original ch sys
- add extra equations for mRNA and conductance regulation 
"""

function (o::OLearyCalcRegulation)(ch::FlowChannel, n::Neuron)
    sys = ch(n)

    (var_conductance,), subst_eqs, _parameters, _defaults, _observed = param_to_state(sys, ch, conductances)

    plast_ch = PlasticisedChannel(ch, o)
    (Ca,) = instantiate_hooks(n, plast_ch, sensed_ions)
    @variables mRNA(t)
    τmRNA, Ca_tgt, τg = get_parameters(o)

    extended_eqs = [D(mRNA) ~ (Ca_tgt - Ca) / τmRNA,
        D(var_conductance) ~ (mRNA - var_conductance) / τg]
    new_eqs = vcat(subst_eqs, extended_eqs)

    new_defaultmap = merge(_defaults, Dict(mRNA => 0.0), default_params(o))

    new_states, new_params = vardivide(Ca, var_conductance, mRNA, τmRNA, Ca_tgt, τg)
    states = [ModelingToolkit.states(sys)..., new_states...]
    parameters = [_parameters..., new_params...]
    plast_sys = ODESystem(new_eqs, t, states, parameters; observed=_observed, defaults=new_defaultmap, name=Symbol(get_name(ch), :_regul))
    return plast_sys
end

function param_to_state(sys, ch, args::Function...)
    _eqs, _states, _parameters, _defaults, _observed = map(
        (equations, states, parameters, ModelingToolkit.get_defaults, observed)) do f
        f(sys)
    end

    (old_conductance,) = map(args) do fun
        getproperty(sys, shorthand_name.(fun(ch))...; namespace=false)
    end
    (var_conductance,) = instantiate_variables(ch, args...)

    # potential bug to fix: only works for a single element of conducts
    _param_indxs = findall(x -> isequal(old_conductance, x), _parameters)

    ## turn conductances into variables
    substitution = Dict(_parameters[_param_indxs]... .=> var_conductance...)
    subst_eqs = [substitute(eq, substitution) for eq in _eqs]

    deleteat!(_parameters, _param_indxs)
    _defaults[var_conductance] = pop!(_defaults, old_conductance)

    return (var_conductance,), subst_eqs, _parameters, _defaults, _observed
end

get_parameters(::OLearyCalcRegulation) = @parameters τmRNA, Ca_tgt, τg
default_params(o::OLearyCalcRegulation) = Dict(get_parameters(o) .=> (o.τmRNA, o.Ca_tgt, o.τg))
