function get_from(d::Dict{DataType,SpeciesDynamics}, func)
    muddled = d |> values .|> func |> x -> filter(y -> !isnothing(y), x)
    isempty(muddled) && return Vector{Num}[]
    return reduce(vcat, muddled)::Vector{Num}
end

function vardivide(v::Num...)
    states = [filter(!ModelingToolkit.isparameter, v)...]
    params = [filter(ModelingToolkit.isparameter, v)...]
    return states, params
end
