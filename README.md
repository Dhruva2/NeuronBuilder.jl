# NeuronBuilder



How to do:

build_neuron function:
Initialise V, Ca

Voltage channels
Build ODESystem for each channel. observed = current.

neuron.V ~ channel.V for channel in Voltage Channels
neuron.resting_potential = channel.resting_potential (how?)


D(V) = sum(channel.observed) for channel in Voltage Channels, Calcium Channels
D(V) += sum(input for input in pins) (connect pins after building)
D(Ca) = sum(channel.observed[2]) for channel in Calcium Channels


### inputs
create subsystem called inputs
each pin is a subsystem of inputs (ie a subsubsystem)


###################
Make equation only in V and Ca. 

Make subsystem for each channel.
Make subsystem for each input 




eqs =   V ~ sum(channels.Iinp) + sum(inputs.Isyn)
        Ca ~ (1/τCa)*(Ca∞ - Ca)

Ca∞ = 0.05 + 0.94*(ICaS + ICaT)

ICaT = ḡCaT*mCaT^3*hCaT*(ECa(Ca) - V)
ICaS = ḡCaS*mCaS^3*hCaS*(ECa(Ca) - V)


### calcium currents
their channel dynamics depend on calcium and voltage, so have to pin both to neuron
their channel dyamics affect calcium and voltage, so same
issue: can't separate the effect of their dynamics on calcium out? or can we? yes we can

NB: makke 2 observed quantities for calcium currents: their contributions to D(V) and D(Ca) respectively

tauCa is just a parameter for the neuron itself

So...
contribution of caS is 1/2τCa)*(Ca∞s - Ca) 
where Ca∞s = 0.05 + 0.94 ICaS

and similar for CaT