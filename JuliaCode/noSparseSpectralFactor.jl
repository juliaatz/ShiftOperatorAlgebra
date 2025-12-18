using Combinatorics
include("ShiftOperator.jl")

function containsTV(A :: Matrix{ShiftOperator})
    n,m = size(A)
    bool = false
    for i = 1:n
        for j = 1:m
            aij = A[i,j]
            if isTV(aij)
                bool = true
                return bool
            end
        end
    end
    return bool
end


function isTV(op)
    coeff = getCoefficients(op)
    coeffTV = coeff[2:end, 2:end]
    if sum(abs.(coeffTV)) < 10^(-10)
        bool = false
    else
        bool = true
    end
    return bool
end
# Set up operators
q = ShiftOperator("q")
id = ShiftOperator(1)
z = ShiftOperator(0)
r = 0.5;                # Discount factor, the value does not affect the result

# Set up matrix
M_base = [-id z z z;
          q'*r -id z z;
          z q'*r q'*r z;
          z z -id q'*r;
          z z z -id]

# plotGraph(M_base, -id, q'*r)  #Uncomment to plot graph

n,m = size(M_base)

# Find permutations
perms = collect(permutations(1:m))
nbr_permutations = length(perms)

nbr_sparseTimeInvariant = 0
nbr_sparseTimeVarying = 0
nbr_nonsparseTimeInvariant = 0
nbr_nonsparseTimeVarying = 0
nbr_runningCholesky = 0

for p = 1:nbr_permutations
    # Set up matrix permutation
    permutation = perms[p]
    M = M_base[:, permutation]
    N = M'*M

    # Find Cholesky
    try
        L = chol(M)

        # Double check factorisation
        E = L*L' - M'*M
        if ~onlyZeros(round.(E, 5))
            error("Factorisation is wrong")
        end

        # Check if sparsity was preserved
        isSparse = true
        for i = 1:m
            for j = 1:i
                nij = N[i,j]
                lij = L[i,j]
                if nij == z && lij != z
                    isSparse = false
                end
            end
        end

        # Check if factorisation is time invariant
        isTI = !containsTV(L)

        if isSparse && isTI
            global nbr_sparseTimeInvariant += 1
        elseif isSparse
            global nbr_sparseTimeVarying += 1
        elseif !isSparse && isTI
            global nbr_nonsparseTimeInvariant += 1
        else
            global nbr_nonsparseTimeVarying += 1
        end
        global nbr_runningCholesky += 1
    catch
    end
end
@show nbr_permutations
@show nbr_runningCholesky
@show nbr_sparseTimeInvariant
@show nbr_sparseTimeVarying
