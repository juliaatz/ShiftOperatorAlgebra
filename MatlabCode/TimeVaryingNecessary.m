% Set up shit operators
clear; close all
q = shiftOperator('q');
one = shiftOperator(1);
z = shiftOperator(0);
r = 0.5;                % Discount factor, the value does not affect the result

% Set up matrix
% Matrix representing the transportation routs and delays
M_base = [-one z z z;
    q'*r -one z z;
    z q'*r q'*r z;
    z z -one q'*r;
    z z z -one];

% plotGraph(M);     % Uncomment if you want to plot the graph

[n,m] = size(M_base);

permutations = perms(1:m);
[nbr_permutations,~] = size(permutations);



nbr_sparseTimeInvariant = 0;
nbr_sparseTimeVarying = 0;
nbr_nonsparseTimeInvariant = 0;
nbr_nonsparseTimeVarying = 0;
nbr_runningCholesky = 0;
for p = 1:nbr_permutations
    % Set up matrix permutation
    permutation = permutations(p,:);
    M = M_base(:,permutation);
    N = M'*M;

    % Find Cholesky
    try
        L = chol(M);
        % Check if sparsity was preserved
        isSparse = true;
        for i = 1:m
            for j = 1:i
                nij = N(i,j);
                lij = L(i,j);
                if nij == z && lij ~= z
                    isSparse = false;
                end
            end
        end
        
        % See if factorisation is time invariant
        isTI = ~containsTV(L);
        
        if isSparse && isTI
            nbr_sparseTimeInvariant = nbr_sparseTimeInvariant + 1;
        elseif isSparse
            nbr_sparseTimeVarying = nbr_sparseTimeVarying + 1;
        elseif ~isSparse && isTI
            nbr_nonsparseTimeInvariant = nbr_nonsparseTimeInvariant + 1;
        else
            nbr_nonsparseTimeVarying = nbr_nonsparseTimeVarying + 1;
        end
        % Double check factorisation
        E = L*L' - M'*M;
        if ~onlyZeros(round(E,12))
            error("Factorisation is wrong")
        end
        nbr_runningCholesky = nbr_runningCholesky + 1;
    catch
    end   
    
end

nbr_permutations
nbr_runningCholesky
nbr_sparseTimeInvariant
nbr_sparseTimeVarying


function bool = containsTV(A)
    [n,m] = size(A);
    bool = false;
    for i = 1:n
        for j = 1:m
            aij = A(i,j);
            if isTV(aij)
                bool = true;
                return
            end
        end
    end
end


function bool = isTV(op)
    coeff = getCoefficients(op);
    coeffTV = coeff(2:end, 2:end);
    if sum(abs(coeffTV),'all') < 10^(-10)
        bool = false;
    else
        bool = true;
    end
end
