# NeuronBuilder

This is a quick package for building small networks of detailed, conductance-based neurons out of ion channels and synapses. The idea is that it's easy to use this template and add your own ion channels /synapses , with your choice of dynamics. 



# Installation and Usage

1. Download this repo
2. Switch to the package manager using `]` and then use `activate .` and `update`
3. Then exit the package manager (`ctrl+c`) and type `using NeuronBuilder`
4. Switch to the shell using `;` and then navigate to the scripts folder
5. Run the script you want using `include XXX.jl`

## Running individual neurons

- Run `include("individual_STG_neurons.jl")` which will make your neurons
- Solve some neuron using ` sol = solve(probAB1, Rodas5())`
- You can convert the solution to an array using `A = Array(sol)` and plot it using `Plots.plot(sol.t,A[1,:])`
