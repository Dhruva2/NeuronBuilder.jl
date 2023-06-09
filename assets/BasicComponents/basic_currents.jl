

function basic_mh_current(S::Species, _m::Pair{Symbol,In}, _h::Pair{Symbol,In}; output=Voltage()::Species) where {In<:Integer}
    
    spec = typeof(S)
    actuated = typeof(output)
    return function current(c::Component, vars...)
        I, E, V, g = find_from([Current{actuated}(), Reversal{spec}(), Voltage(), Conductance{spec}()], vars...)

        m,h = map((_m, _h)) do el
            (last(el) > 0) ? (return find_from(UntrackedQuantity(first(el)), vars...)) : (return 1.0)
        end

        # isempty(m) && m = 1.0
        # isempty(h) && h = 1.0 
        return [
            I ~ g * m^last(_m) * h^last(_h) * (E - V),
        ]
    end
end


function multiple_ion_current(sensed::Vector{S}) where {S<:Species}
    
end
