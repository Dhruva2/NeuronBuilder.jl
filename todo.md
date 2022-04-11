1. Consistent way of turning parameters into states with dynamics. Regulated is one example. More generally want it to work to add synapses with dynamics.

Idea: PlasticityRule(Component). Regulated is a plasticity rule on ion channels. EG stdp is a plasticity rule on synapses.

Each plasticity rule just turns parameters into states with specified dynamics

2. ComponentInfo abstract type
subtypes are the functions like ionic_conductance which get info from the components. Then can get them to behave in common ways



3. Find a nice way of adding callbacks. need to access e.g. the index of the voltage of the subsystems of different neurons in the solution array

VoltageIndex(ODESystem, i)  = index of voltage of ith neuron


4. ask andrea to get rid of the 14.96 and the prinz_conversion etc in the package. Just make an overall constant OUTSIDE of the package that multiplies the conductances. Will need different constants for calcium and voltage?

### ROugh
Units:
if e.g. maximal conductances is divided by area, then multiply by internal soma area. otherwise don't


Plasticity Rule:

1. Finds relevant parameters


"""
Ion channel logic:

struct Ion
Channel{Vector{Ion}} # which ions it senses and/or actuates
Reversal(channel) = Reversal{Ions}
All channels sense and actuate electricity


"""

reversals(ECa) can be a function

flows are the ions that are not voltage that are quantified by the neurons
CaS channels:

so flows(sodium channel) = :na 

When the soma sees the channel, it only adds state equations for the flows that it senses

So the LiuNeuron would be Neuron{Calcium}
and then it would add a state equation for calcium and include the flow
the soma itself would add a constant such as 14.96 to normalise the flow based on the ionic units in the soma.


All neurons track voltage
so don't have neuron{Voltage, Calcium}


new build neuron:
instead of summed calcium flux: find intersection of the ionic flows with the ions recorded in the soma.  which is just calcium here.
Then have an external update type of equation for voltage and current in the soma, particular to the type of neuron

ABNeuron = Neuron{Calcium}(StandardVoltageUpdate(), LiuCalciumUpdate(neuron flux)) 
If the requisite arguments are not entered, default to a stnadard update:
1 / tau_ion (or capacitance for voltage) * summed flux

LiuCalciumUpdate gives the RHS of D(Ca), given an input summed_flows for calcium

MOVE ECa to property of soma, not cas and cat channels


## Add convenience function
which takes in a channel and checks it for consistency:
    its ODESystem has a voltage 


which takes in a geometry and checks it for consistency:
    it calculates reversals
    it calculates capacitance