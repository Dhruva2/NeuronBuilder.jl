module Liu

using ModelingToolkit
using ..NeuronBuilder


include("channel_dynamics.jl")
include("channels.jl")
include("species_dynamics.jl")
include("synapses.jl")
end