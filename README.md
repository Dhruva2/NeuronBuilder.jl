# NeuronBuilder

This is a quick package for building small networks of detailed, conductance-based neurons out of ion channels and synapses. The idea is that it's easy to use this template and add your own ion channels /synapses , with your choice of dynamics. 

Issue:

When I connect up the STG neurons, I get an instability as soon as they start to spike. What happens? 

For the synapse internal state variables, |ds/dt| -> infinity.
Why? Because it is proportional to 1/τs, and τs -> 0.
When does τs -> 0. ? When s̄(Vpre) = 1. -> When exp(Vth - Vpre) is large.
