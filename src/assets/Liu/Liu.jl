module Liu

using ModelingToolkit
using ..NeuronBuilder

import ..channel_dynamics, ..ionic_current, ..calcium_current, ..external_params, ..flow_variables

import ..ions
import ..reversals, ..currents, ..conductances

# import ..Flows, ..InstantiateVariables, ..IonicCurrent

include("channels.jl")

export channel_dynamics, ionic_current, ionic_conductance, calcium_current, external_params

export CalciumFlow, IonicConductance, IonicCurrent
end

"""
Deliberately not exporting ion channel names so they don't conflict.

eg Na can only be accessed as Liu.Na from NeuronBuilder

Deliberately exporting other structs that act on components. because it is one struct, not defined separaetly in different modules, that acts on module-specific structs
"""