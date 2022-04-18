module Liu

using ModelingToolkit
using ..NeuronBuilder

import ..get_parameters, ..get_states, ..default_params, ..default_states

include("channels.jl")
include("calc_dynamics.jl")

export get_parameters, get_states, default_params, default_states

end

"""
Deliberately not exporting ion channel names so they don't conflict.

eg Na can only be accessed as Liu.Na from NeuronBuilder

Deliberately exporting other structs that act on components. because it is one struct, not defined separaetly in different modules, that acts on module-specific structs
"""