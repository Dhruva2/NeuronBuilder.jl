# NeuronBuilder

This is a quick package for building small networks of detailed, conductance-based neurons out of ion channels and synapses. The idea is that it's easy to use this template and add your own ion channels /synapses, with your choice of dynamics. 

If you want a more flexible platform to build neuron models from basic components you should check out the more comprehensive package [Conductor.jl](https://github.com/wsphillips/Conductor.jl). 


# Installation and Usage

1. Download this repo
2. Switch to the package manager using `]` and then use `activate .` and `update`
3. Then exit the package manager (`ctrl+c`) and type `using NeuronBuilder`
4. Switch to the shell using `;` and then navigate to the scripts folder
5. Run the script you want using `include XXX.jl`

## Running individual neurons

- Build channels with `Liu.Na(g)` or `Prinz.Na(g)` if your sodium conductance has value `g`. 
- To reproduce sets of channels as reported in the original papers (https://www.jneurosci.org/content/18/7/2309) (https://journals.physiology.org/doi/full/10.1152/jn.00641.2003) use the conversion factors for your choice of module: Liu to Neuron Builder -> `Liu_conversion` or Prinz to NeuronBuilder -> `Prinz_conversion`. These are specifically given in the individual neuron script.
- From the Julia REPL run `include("individual_neurons.jl")` which will make your neurons

## Running a network of neurons

- The `connected_STG.jl` script shows how to add synapses between neurons.
- The group of neurons and their connections is an `ODESystem`; the names of the variables in this system can be listed with `@show final_sys.states`. The problem that takes this system is a regular `ODEProblem` and You can convert the solution to an array using `A = Array(sol)` and plot it using `Plots.plot(sol.t,A[1,:])`. 


