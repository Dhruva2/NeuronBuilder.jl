function build_group(neurons::Vector{N}; name=:group) where {N<:Any}
    eqs = Array{Equation,1}(undef, 0)
    ODESystem(eqs, t, [], []; name=name, systems=[el.sys for el in neurons])
end

function add_connection(group, pre, post, syn::EmptyConnection; kwargs...)
    return group
end

function add_connection(group, pre_n, post_n, syn::Synapse; name=ModelingToolkit.getname(group), i=1)
    pre = pre_n.sys
    post = post_n.sys
    prename, postname = ModelingToolkit.getname.([pre, post])
    synapse_sys = syn(; name=Symbol(prename, :to, postname, get_name(syn)))

    oldeqs = ModelingToolkit.get_eqs(group)
    neweqs = [
        get_pre(syn, pre) ~ my_pre(synapse_sys),
        get_post(syn, post) ~ my_post(synapse_sys),
        syn_current(syn, synapse_sys.sys) ~ getproperty(post, Symbol(post_connector(syn), i))
    ]

    eqs = cat(oldeqs, neweqs, dims=1)
    systems = cat(ModelingToolkit.get_systems(group), synapse_sys.sys, dims=1)

    return connected = ODESystem(eqs, t, [], [];
        name=name,
        systems=systems)
end

function in_degree(edge_fn, num_nodes, j)
    sum(1:num_nodes) do i
        length(edge_fn(i, j)) - (typeof(edge_fn(i, j)) == Tuple{EmptyConnection})
    end
end

function build_network(nodefn, edgefn, which_nodes::AbstractRange; name=:network)
    eqs = Array{Equation,1}(undef, 0)
    num_nodes = length(which_nodes)
    neurons = [nodefn(j)(in_degree(edgefn, num_nodes, j)) for j in 1:num_nodes]
    network = ODESystem(eqs, t, [], []; name=name, systems=[el.sys for el in neurons])
    fill_counter = ones(Int64, num_nodes)
    for i in which_nodes
        for j in which_nodes
            for syn in edgefn(i, j)
                network = add_connection(network, neurons[i], neurons[j], syn; i=fill_counter[j])
                !(typeof(syn) == EmptyConnection) && (fill_counter[j] += 1)
            end
        end
    end
    return network |> structural_simplify
end
