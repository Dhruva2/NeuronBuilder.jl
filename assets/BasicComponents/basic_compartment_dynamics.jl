"""
i.e. C dv/dt = sum(currents{species})

- where C = Capacitance{species} (species = voltage usually, but can also be pseudo-capacitance for ions)

"""
function basic_influx_dynamics(S::Species)
    function f(c::Component, owned_systems::Vector{S}, internal_vars...) where S <: ModelingToolkit.AbstractSystem
        V, C = find_from([Voltage(), Capacitance{Voltage}()], internal_vars...)
        flux = sum(owned_systems) do el
            (hasproperty(el, shorthand_name(Current{Voltage}()))) ? (getproperty(el, shorthand_name(Current{Voltage}()))) : 0.0
        end
        return D(V) ~ (1 / C) * (flux)
    end
end

