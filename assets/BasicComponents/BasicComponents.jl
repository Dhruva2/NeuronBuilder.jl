module BasicComponents

using ModelingToolkit
using ..NBuilder

import ..sensed, ..actuated, ..is_dynamic, ..has_dynamics, ..tagged_internal_variables, ..untagged_internal_variables


include("basic_currents.jl")
include("basic_ion_channels.jl")
include("basic_compartment_dynamics.jl")

export sensed, actuated, is_dynamic, has_dynamics, tagged_internal_variables, untagged_internal_variables
end