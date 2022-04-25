module HH

using ModelingToolkit
using ..NeuronBuilder

import ..get_parameters, ..get_states, ..default_params, ..default_states

include("channels.jl")

export get_parameters, get_states, default_params, default_states

end