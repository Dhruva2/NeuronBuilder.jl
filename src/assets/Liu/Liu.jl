module Liu

using ModelingToolkit
#include("../../channel_types.jl")

using ...NeuronBuilder
include("channels.jl")

export channel_dynamics


end