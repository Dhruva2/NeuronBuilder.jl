module Prinz

using ModelingToolkit
using ..NeuronBuilder

import ..channel_dynamics
include("channels.jl")

export channel_dynamics

end