# NeuronBuilder

This is a quick package for building small networks of detailed, conductance-based neurons out of ion channels and synapses. The idea is that it's easy to use this template and add your own ion channels / synapses, with your choice of dynamics. Iteratively adding these components to a neuron is done using ModelingToolkit, which makes it scalable to any number of ion channels and synapses.

If you want a more flexible platform to build neuron models, with e.g. multiple compartments, from basic components you should check out the more comprehensive package [Conductor.jl](https://github.com/wsphillips/Conductor.jl) (in active development). 


# Installation and Usage

NeuronBuilder is available as a registered package and has a tagged release (v0.1.0)
```
#From Julia REPL
] add NeuronBuilder
```
Once you exit the package manager (`ctrl+c`), type `using NeuronBuilder`.

To try out the demo scripts:
1. From a terminal, clone this repository
```
git clone https://github.com/Dhruva2/NeuronBuilder.jl
cd NeuronBuilder.jl
git checkout v0.2.4
```
2. Navigate to the scripts folder. Once there, activate the environment 
```
] activate .
```
3. Open a Julia session and run the script you want, for example `include("neuron_liu.jl")`

## Running individual neurons

- A sodium channel is created with `Liu.Na(g)` or `Prinz.Na(g)` if the conductance has value `g`. 
- To reproduce sets of channels as reported in the papers [Liu et. al. 1998](https://www.jneurosci.org/content/18/7/2309) and [Prinz et. al. 2003](https://journals.physiology.org/doi/full/10.1152/jn.00641.2003) use the conversion factors `Liu_conversion` and `Prinz_conversion` to get the right units. You can see how these are specifically given in the scripts and are multiplying with the original `g` value from the papers. 
- The `g_values.jl` file has a small collection of conductances coming from various sources. You can copy-paste any of these into the script that simulates a single neuron.

## Running a network of neurons

- The `connected_STG.jl` script shows how to add synapses between neurons and reproduces the triphasic rhythm of the STG found in [Prinz et. al. 2004](http://www.nature.com/articles/nn1352).
- Synapses also get a conversion factor which depends on the geometry of the somatic compartment.

## Adding your own libraries of custom ion channels

- Fork NeuronBuilder 
- Go to the src/assets folder, you'll find two modules called Liu and Prinz. 
- Each module has a list of channels with specified dynamics
- You can copy-paste the structure of those modules, just changing the dynamics and channel names.
- Feel free to send us the library you built to push to the master branch! :)

## Acknowledgements
This work was funded by European Research Council grant 716643 FLEXNEURO and HFSP grant RGY0069/2017, (Principal Investigator Timothy O’Leary). 
