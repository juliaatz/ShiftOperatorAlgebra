function M2ss(M)
    q = ShiftOperator("q")
    one = ShiftOperator(1)

    p,m = size(M)

    # Ensure that M has the correct structure
    M0 = monomialCoefficients(M,1,1)
    M1 = monomialCoefficients(M,2,1)
    if M0*one + M1*q' != M
        throw(DomainError(M,"the matrix M does not fullfil the requirments"))
    end

    n = p + m

    A = zeros(n,n)
    B = zeros(n,m)
    C = zeros(p,n)
    D = zeros(p,m)

    # Add information about storage and what the output states are
    oddStates = 1:2:n
    for i = oddStates
        A[i,i] = 1
        C[Int8((i+1)/2), i] = 1
    end

    # Add information about the flow
    for i = 1:m
        edge = M[:,i]
        # Find target and source node
        connections = findNonZeroElements(edge)
        connection1 = connections[1]
        connection2 = connections[2]

        if edge[connection1] == -one
            sourceNode = connection1
            targetNode = connection2
        elseif edge[connection2] == -one
            sourceNode = connection2
            targetNode = connection1
        else
            throw(DomainError(M, "the matrix M does not fulfill the requiements"))
        end

        # Delayed state
        delayedState = 2*i
        B[delayedState, i] = 1

        # Delayed inflow
        targetState = 2*targetNode - 1
        A[targetState, delayedState] = 1

        # Outflow
        sourceState = 2*sourceNode - 1
        B[sourceState, i] = -1
    end

    return A,B,C,D
end

function findNonZeroElements(V)
    n = length(V)

    z = ShiftOperator(0)
    is = []
    for i = 1:n
        vi = V[i]
        if vi != z
            push!(is,i)
        end
    end
    return is
end


