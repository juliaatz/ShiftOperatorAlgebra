function plotGraph(M, outflow_pattern, inflow_pattern)
    A, edgelabel_dict = adjecencyM(M, inflow_pattern, outflow_pattern)
    G = DiGraph(A)

    p,m = size(M)
    graphplot(G, names=1:p, edgelabel=edgelabel_dict, curves=false, nodeshape=:circle)
end

function adjecencyM(M, inflow_pattern, outflow_pattern)
    p,m = size(M)
    
    A = Base.zeros(p,p)
    edgelabel_dict = Dict()
    for i = 1:m
        inflow = -1
        outflow = -1
        for j = 1:p
            op = M[j,i]
            if op == outflow_pattern
                outflow = j
            end

            if op == inflow_pattern
                inflow = j
            end
        end
        try
            A[outflow, inflow] = 1
            edgelabel_dict[(inflow,outflow)] = string(i)
        catch
            throw(DomainError(M, "M does not have one inflow and one outflow per edge"))
        end
    end

    return A, edgelabel_dict
end


