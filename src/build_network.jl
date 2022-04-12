"""
1. number of hooks for each neuron. This is in-degree of neuron
"""




function build_network(nodefn, edgefn, num_nodes; name=:network)
    eqs = Array{Equation,1}(undef, 0)
    neurons = [nodefn(j)(in_degree(edgefn, num_nodes, j)) for j in 1:num_nodes]
    network = ODESystem(eqs, t, [], []; name=name, systems=[el.sys for el in neurons])
    fill_counter = ones(Int64, num_nodes)
    for i = 1:num_nodes
        for j = 1:num_nodes
            for syn in edgefn(i, j)
                network = add_connection(network, neurons[i], neurons[j], syn; i=fill_counter[j])
                !(typeof(syn) == EmptyConnection) && (fill_counter[j] += 1)
            end
        end
    end
    return network |> structural_simplify
end


function in_degree(edge_fn, num_nodes, j)
    sum(1:num_nodes) do i
        length(edge_fn(i, j)) - (typeof(edge_fn(i, j)) == Tuple{EmptyConnection})
    end
end
# if typeof(syn) !== EmptyConnection
#     println(typeof(syn), i, j)
#     println(ModelingToolkit.getname(neurons[i].ODESystem), " is pre and ", ModelingToolkit.getname(neurons[j].ODESystem), " is post \n")
# end
