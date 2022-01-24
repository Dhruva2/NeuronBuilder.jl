######### helper functions ###########

get_g(ch::IonChannel) = fieldnames(typeof(ch))[1]
get_R(ch::RegIonChannel) = fieldnames(typeof(ch))[2]

"""
return unparameterised type as symbol
"""
function get_name(ch::Component)
    Base.typename(ch |> typeof).name |> Symbol
end

calcium_current(c::IonChannel, s::ODESystem) = Num(0)
voltage_hook(V, cs) = V ~ cs.sys.V
calcium_hook(Ca, cs) = Ca ~ cs.sys.Ca

isLeak(c::IonChannel) = typeof(c) in [NeuronBuilder.Liu.Leak{Float64}, NeuronBuilder.Prinz.Leak{Float64}]
# """
# does the channel have a reversal
# """
# get_reversal(c::Conductance) = nothing