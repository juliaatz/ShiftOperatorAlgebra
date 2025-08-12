function [fig, xPos, yPos] = plotGraph(M)
%PLOTGRAPH(M) plots the graph corresponding to the matrix of shift
%operators M


% Build adjecency matrix
[n, m] = size(M);

startNodes = zeros(1,m);
endNodes = zeros(1,m);
edgeLabels = 1:m;

for i = 1:m
    
    for j = 1:n
        element = M(j,i);
        coeff = getCoefficients(element);
        % Find outflow
        direktFlow = coeff(1,1);
        if( direktFlow < 0 )
            outFlow = j;
        end
        
        % find inflow
        try
            delayedFlow = coeff(2,1);
            if ( delayedFlow > 0 )
                inFlow = j;
            end
        catch
        end
    end

    startNodes(i) = outFlow;
    endNodes(i) = inFlow;
end
graph = digraph(startNodes, endNodes, edgeLabels);

fig = figure;
plotData = plot(graph, 'Layout', 'force', 'EdgeLabel', graph.Edges.Weight, 'ArrowSize', 17, 'LineWidth', 2, 'NodeColor', 'k', 'MarkerSize', 7);
xPos = plotData.XData;
yPos = plotData.YData;
end

