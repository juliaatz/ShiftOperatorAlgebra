function M = randomTree(m, mixedDirection, outflow_pattern, inflow_pattern)
%RANDOMTREE creates a matrix of shift operators corresponding to a random
%tree graph
%
%RANDOMTREE(m, mixedDirection, outflow_pattern, inflow_pattern)
%   m - The number of edges in the graph
%   mixedDirection - If true (false) the direction of the edges will be
%                       random (deterministic)
%   outflow_pattern - The element in the matrix to represent outflow
%   inflow_pattern - The element in the matrix to represent inflow
    p = m + 1;
    M = shiftOperator.zeros(p,m);
    M_unidirect = shiftOperator.zeros(p,m);
    for i = 1:m
        % Add edge
        base_node = randi([1,i],1); % Random base node
        new_node = i + 1;
        
        random_direction = randi([0,1],1);
        if ~mixedDirection
            random_direction = 1;
        end
        
        if random_direction
            M(base_node, i) = outflow_pattern;
            M(new_node, i) = inflow_pattern;
        else
            M(new_node, i) = outflow_pattern;
            M(base_node, i) = inflow_pattern;
        end
        
        if mixedDirection
            M_unidirect(base_node, i) = outflow_pattern;
            M_unidirect(new_node, i) = inflow_pattern;
        end
    end
    M = M(end:-1:1,end:-1:1);
    M_unidirect = M_unidirect(end:-1:1,end:-1:1);

    % Reorder edges and nodes to make it more logical
    source_node = p;
    if mixedDirection
        [improved_rev_orderV, improved_rev_orderE] = findRevOrder(M_unidirect, source_node, outflow_pattern, inflow_pattern);
    else
        [improved_rev_orderV, improved_rev_orderE] = findRevOrder(M, source_node, outflow_pattern, inflow_pattern);
    end
    improved_orderE = improved_rev_orderE(end:-1:1);
    improved_orderV = [improved_rev_orderV(end:-1:1), p];
    M = M(:, improved_orderE);
    M = M(improved_orderV,:);
end

function [orderV, orderE] = findRevOrder(M, base_node, outflow_pattern, inflow_pattern)
    % Find connected edges
    [~, edges] = findPattern(M(base_node,:), outflow_pattern);

    % Add edges
    orderE = edges;
    orderV = [];
    % Add nodes
    for e = edges
        % Find new base node
        outflow_node = findPattern(M(:,e), inflow_pattern);
        orderV = [orderV, outflow_node];
    end

    % Recursive part
    for e = edges
        % Find new base node
        outflow_node = findPattern(M(:,e), inflow_pattern);
        
        % find order from there
        [subgraph_orderV, subgraph_orderE] = findRevOrder(M, outflow_node, outflow_pattern, inflow_pattern);
        
        % Append order
        orderE = [orderE, subgraph_orderE];
        orderV = [orderV, subgraph_orderV];
    end
end

function [is,js] = findPattern(M, op)
    % returns all entry indecies in M that are equal to op
    [n,m] = size(M);
    is = [];
    js = [];
    for i = 1:n
        for j = 1:m
            mij = M(i,j);
            if mij == op
                is(end+1) = i;
                js(end+1) = j;
            end
        end
    end
end
