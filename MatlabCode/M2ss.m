function sys = M2ss(M)
    %M2SS Finds a state space representation from the a matrix of shift operators.
    %
    %In M2SS(M), the edges in M are assumed to contain the elements
    %(-1, r*q') and the remaining dynamics in the system are assumed to be
    %captured by the dynamics are 1/(1+r*q').
    %   The even states correspond to control signals x_2i = u_i, and the
    %odd states to output signals x_2i+1 = y_i-1 

    q = shiftOperator('q');
    one = shiftOperator(1);
    
    [p,m] = size(M);
    
    % Ensure that there are no longer delays than one-step delays
    M0 = monomialCoefficients(M,1,1);
    M1 = monomialCoefficients(M,2,1);
    if (M0 + M1*q'~= M)
        errorStruct.message = "The matrix M does not fullfil the requirments";
        errorStruct.identifier = "M2ss:Mstructure";
        error(errorStruct);
    end
 
    n = p + m;
    
    A = zeros(n);
    B = zeros(n,m);
    C = zeros(p,n);
    D = zeros(p,m);
    
    % Add information about storage and what are the output states
    odd_states = 1:2:n;
    for i = odd_states
        A(i,i) = 1;
        C((i+1)/2, i) = 1;
    end 

    % Add information about the flow
    for i  = 1:m
        edge = M(:,i);
        % Find target and source node
        connections = find(edge);
        connection1 = connections(1);
        connection2 = connections(2);
        
        if edge(connection1)== -one
            sourceNode = connection1;
            targetNode = connection2;
        elseif edge(connection2) == -one
            sourceNode = connection2;
            targetNode = connection1;
        else
            errorStruct.message = "The matrix M does not fullfil the requirments";
            errorStruct.identifier = "M2ss:Moutflow";
            error(errorStruct)
        end
        
        % Delayed state
        delayState = 2*i;
        B(delayState, i) = 1;
        
        % Delayed inflow
        targetState = 2*targetNode - 1;
        A(targetState, delayState) = 1;
        
        % Ouflow
        sourceState = 2*sourceNode - 1;
        B(sourceState, i) = - 1;
        
    end

    sys = ss(A,B,C,D);
end