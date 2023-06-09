function calcium_dynamics(; τCa=20.0, Ca∞=0.05, calc_multiplier=14.96 * 0.0628)
    return function f(c::Component, owned_systems::Vector{S}, internal_vars...) where {S<:ModelingToolkit.AbstractSystem}
        Ca, Cₘ = find_from([Calcium(), Capacitance{Voltage}()], internal_vars...)
        flux = sum(owned_systems) do el
            (hasproperty(el, shorthand_name(Current{Calcium}()))) ? (getproperty(el, shorthand_name(Current{Calcium}()))) : 0.0
        end
        return D(Ca) ~ (1 / τCa) * (-Ca + Ca∞ + (calc_multiplier * flux / Cₘ))
    end
end

function calcium_reversal_dynamics()
    function (c::Component, owned_systems::Vector{S}, internal_vars...) where {S<:ModelingToolkit.AbstractSystem}
        Ca, ECa = find_from([Calcium(), Reversal{Calcium}()], internal_vars...)
        return ECa ~ (500.0) * (8.6174e-5) * (283.15) * (log(max((3000.0 / Ca), 0.001)))
    end
end