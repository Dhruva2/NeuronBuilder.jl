module Liu

using ModelingToolkit
#include("../../channel_types.jl")

using ..NeuronBuilder
import ..channel_dynamics
include("channels.jl")

export channel_dynamics

end