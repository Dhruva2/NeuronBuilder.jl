module Liu

using ModelingToolkit
using ..NBuilder


include("channel_dynamics.jl")
include("channels.jl")
include("species_dynamics.jl")
end