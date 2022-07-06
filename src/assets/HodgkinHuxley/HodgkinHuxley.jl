module HodgkinHuxley

using ModelingToolkit
using ..NeuronBuilder

import ..channel_dynamics
import ..get_parameters, ..get_states, ..default_params, ..default_states

include("channels.jl")
include("synapses.jl")

export channel_dynamics
export get_parameters, get_states, default_params, default_states

end
