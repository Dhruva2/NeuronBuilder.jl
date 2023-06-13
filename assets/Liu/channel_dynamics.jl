
function leaky_integrator(V, g, τ, g∞)
    return D(g) ~ (1 / τ(V)) * (g∞(V) - g)
end

gating_dynamics = Dict(
    :Na => Dict(
        :m∞ =>  V -> 1.0 / (1.0 + exp((V + 25.5) / -5.29)),
        :τm =>  V  -> 1.32 - 1.26 / (1 + exp((V + 120.0) / -25.0)),
        :h∞ => V -> 1.0 / (1.0 + exp((V + 48.9) / 5.18)),
        :τh =>   V -> (0.67 / (1.0 + exp((V + 62.9) / -10.0))) * (1.5 + 1.0 / (1.0 + exp((V + 34.9) / 3.6)))
    ),
    :CaS => Dict(
        :m∞ => V -> 1.0 / (1.0 + exp((V + 33.0) / -8.1)),
        :τm => V -> 1.4 + 7.0 / (exp((V + 27.0) / 10.0) + exp((V + 70.0) / -13.0)),
        :τh => V ->  60.0 + 150.0 / (exp((V + 55.0) / 9.0) + exp((V + 65.0) / -16.0)),
        :h∞ => V -> 1.0 / (1.0 + exp((V + 60.0) / 6.2))
    ),
    :CaT => Dict(
        :m∞ => V -> 1.0/(1.0+exp((V + 27.1) / -7.2)),
        :τm => V -> 21.7 - 21.3 / (1.0 + exp((V + 68.1) / -20.5)),
        :τh => V -> 105.0 - 89.8 / (1.0 + exp((V + 55.0) / -16.9)),
        :h∞ => V -> 1.0 / (1.0 + exp((V + 32.1) / 5.5))
    ),
    :KCa => Dict(
        :m∞ => (V, Ca) -> (Ca / (Ca + 3.0)) / (1.0 + exp((V + 28.3) / -12.6)),
        :τm => V -> 90.3 - 75.1 / (1.0 + exp((V + 46.0) / -22.7))
    ),
    :Ka => Dict(
        :m∞ => V -> 1.0 / (1.0 + exp((V + 27.2) / -8.7)),
        :τm => V -> 11.6 - 10.4 / (1.0 + exp((V + 32.9) / -15.2)),
        :τh => V -> 38.6 - 29.2 / (1.0 + exp((V + 38.9) / -26.5)),
        :h∞ => V -> 1.0 / (1.0 + exp((V + 56.9) / 4.9))
    ),
    :Kdr => Dict(
        :m∞ => V -> 1.0 / (1.0 + exp((V + 12.3) / -11.8)),
        :τm => V -> 7.2 - 6.4 / (1.0 + exp((V + 28.3) / -19.2))
    ),
    :H => Dict(
        :m∞ => V -> 1.0 / (1.0 + exp((V + 70.0) / 6.0)),
        :τm => V -> (272.0 + 1499.0 / (1.0 + exp((V + 42.2) / -8.73)))
    )
)



synaptic_gating_dynamics = Dict(
    :Chol => Dict{Symbol,Any}(
        :s̄ => (Vth, Vpre, δ) -> 1.0 / (1.0 + exp((Vth - Vpre) / δ))
    ),
    :Glut => Dict(
        :s̄ => x -> x)
)

push!(synaptic_gating_dynamics[:Chol], :τs => (Vth, Vpre, δ, k₋) -> (1.0 - (synaptic_gating_dynamics[:Chol][:s̄](Vth, Vpre, δ)) / k₋))


channel_dynamics = Dict(
    :Na => Dict(
        :m =>   (c::Component, vars...) ->
                    leaky_integrator(find_from([Voltage(), UntrackedQuantity(:m)], vars...)..., gating_dynamics[:Na][:τm], gating_dynamics[:Na][:m∞])
                ,
        :h =>   (c::Component, vars...) ->
                    leaky_integrator(find_from([Voltage(), UntrackedQuantity(:h)], vars...)..., gating_dynamics[:Na][:τh], gating_dynamics[:Na][:h∞])
    ),
    :CaS => Dict(
        :m => (c::Component, vars...) ->
                    leaky_integrator(find_from([Voltage(), UntrackedQuantity(:m)], vars...)..., gating_dynamics[:CaS][:τm], gating_dynamics[:CaS][:m∞]),
        :h => (c::Component, vars...) ->
            leaky_integrator(find_from([Voltage(), UntrackedQuantity(:h)], vars...)..., gating_dynamics[:CaS][:τh], gating_dynamics[:CaS][:h∞])
    ),
    :CaT => Dict(
        :m => (c::Component, vars...) ->
                    leaky_integrator(find_from([Voltage(), UntrackedQuantity(:m)], vars...)..., gating_dynamics[:CaT][:τm], gating_dynamics[:CaT][:m∞]),
        :h => (c::Component, vars...) ->
            leaky_integrator(find_from([Voltage(), UntrackedQuantity(:h)], vars...)..., gating_dynamics[:CaT][:τh], gating_dynamics[:CaT][:h∞])
    ),
    :KCa => Dict(
        :m =>  (c::Component, vars...) -> begin
                τm, m∞ = gating_dynamics[:KCa][:τm], gating_dynamics[:KCa][:m∞]
                V, Ca, m = find_from([Voltage(), Calcium(), UntrackedQuantity(:m)], vars...)
                D(m) ~ (1 / τm(V)) * (m∞(V, Ca) - m)
            end
    ),
    :Ka => Dict(
        :m => (c::Component, vars...) ->
                leaky_integrator(find_from([Voltage(), UntrackedQuantity(:m)], vars...)..., gating_dynamics[:Ka][:τm], gating_dynamics[:Ka][:m∞]),
        :h => (c::Component, vars...) ->
            leaky_integrator(find_from([Voltage(), UntrackedQuantity(:h)], vars...)..., gating_dynamics[:Ka][:τh], gating_dynamics[:Ka][:h∞])
    ),
    :Kdr => Dict(
        :m => (c::Component, vars...) ->
            leaky_integrator(find_from([Voltage(), UntrackedQuantity(:m)], vars...)..., gating_dynamics[:Kdr][:τm], gating_dynamics[:Kdr][:m∞])
    ),
    :H => Dict(
        :m => (c::Component, vars...) ->
            leaky_integrator(find_from([Voltage(), UntrackedQuantity(:m)], vars...)..., gating_dynamics[:H][:τm], gating_dynamics[:H][:m∞]) 
    )
)

synaptic_channel_dynamics = Dict(
    :s => (c::Component, vars...) -> begin
        g = synaptic_gating_dynamics[:Chol]
        τs = g[:τs]
        s̄ = g[:s̄]
        s, k₋, Vth, δ, E, Vpre, Vpost = find_from([UntrackedQuantity(:s), UntrackedQuantity(:k₋), UntrackedQuantity(:Vth), UntrackedQuantity(:δ), Reversal{Choline}(), Previous{Voltage}(), Voltage()], vars...)
        D(s) ~ (1 / τs(Vth, Vpre, δ, k₋)) * (s̄(Vth, Vpre, δ) - s)
    end
)