module Prinz

using ModelingToolkit
using ..NeuronBuilder

import ..get_parameters, ..get_states, ..default_params, ..default_states

include("channels.jl")
include("calc_dynamics.jl")

export get_parameters, get_states, default_params, default_states

end

