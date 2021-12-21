module Liu

using ModelingToolkit
using ..NeuronBuilder

import ..channel_dynamics, ..ionic_current, ..calcium_current, ..external_params

include("channels.jl")

export channel_dynamics, ionic_current, calcium_current, external_params

end
