function randomTree(m, mixedDirection, outflow_pattern, inflow_pattern)
    p = m+1
    M = shiftopZeros(ShiftOperator,p,m)
    M_unidirected = shiftopZeros(ShiftOperator,p,m)
    for i = 1:m
        # Add edge
        base_node = rand(1:i)
        new_node = i + 1

        random_direction = rand(0:1)

        if !mixedDirection
            random_direction = 1
        end

        if random_direction == 1
            M[base_node,i] = outflow_pattern
            M[new_node,i] = inflow_pattern
        else
            M[new_node,i] = outflow_pattern
            M[base_node,i] = inflow_pattern
        end

        if mixedDirection
            M_unidirected[base_node,i] = outflow_pattern
            M_unidirected[new_node,i] = inflow_pattern
        end
    end
    M = M[end:-1:1, end:-1:1]
    M_unidirected = M_unidirected[end:-1:1, end:-1,1]

    # Reorder edges and nodes to make it more logical
    source_node = p
    if mixedDirection
        improved_rev_orderN, improved_rev_orderE = findRevOrder(M_unidirected, source_node, outflow_pattern, inflow_pattern)
    else
        improved_rev_orderN, improved_rev_orderE = findRevOrder(M, source_node, outflow_pattern, inflow_pattern)
    end
    improved_orderN = improved_rev_orderN[end:-1:1]
    push!(improved_orderN, p)
    improved_orderE = improved_rev_orderE[end:-1:1]
    M = M[:, improved_orderE]
    M = M[improved_orderN, :]
    return M
end

function findRevOrder(M, base_node, outflow_pattern, inflow_pattern)
    # Find connected edges
    edges = findPattern(M[base_node, :], outflow_pattern)

    # Add edges
    orderE = Int[]
    orderN = Int[]
    for e = edges
        # Find new base node
        outflow_node = findPattern(M[:,e], inflow_pattern)
        push!(orderN, outflow_node[1])
        push!(orderE, e)
    end

    # Recursive part
    for e = edges
        # Find new base node
        outflow_node = findPattern(M[:,e], inflow_pattern)
         
        # Find order from new base node
        subgraph_orderN, subgraph_orderE = findRevOrder(M, outflow_node, outflow_pattern, inflow_pattern)

        # Append order
        for ind = subgraph_orderN
            push!(orderN, ind)
        end
        for ind = subgraph_orderE
            push!(orderE, ind)
        end
    end
    return orderN, orderE
end

function findPattern(V, op)
    # Returns all entry indices in M that are equal to op
    n = length(V)
    is = Int[]
    for i = 1:n
        vi = V[i]
        if vi == op
            push!(is, i)
        end
    end
    return is
end
