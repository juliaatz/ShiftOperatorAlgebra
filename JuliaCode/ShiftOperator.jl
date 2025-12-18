include("SignalL2.jl")
struct ShiftOperator
    shiftMatrix     # matrix M for a operate op such that op = [1 q q^2...]M[1 q' (q')^2 ...]'.

    # SHIFTOPERATOR
    # This is a class of shift operators and matrices of shift operator on
    # l2+. These operators are built from the forwardshift operator q and
    # its adjoint the backward shift operator q'.
    #
    # Given a signal y in l2+
    #       q : [y1, y2, y3, ...] -> [y2, y3, y4, ...]
    #       q': [y1, y2, y3, ...] -> [0, y1, y2, ...].
    #
    # This class allows you to perform algebraic operations with these
    # obejcts.

    function ShiftOperator(x)
        # SHIFTOPERATOR Constructs a shift operator from its shift coefficients or a string.
        # SHIFTOPERATOR(shiftData) creates a shift operator.
        # ShiftData can either be a matrix or a string.
        #   
        # Matrix: shiftData is the matrix M of coefficeients such that
        # the operator is [1 q q^2...]M[1 q' (q')^2 ...]'.
        #   
        # String: ShiftData is on of the following 4 strings:
        #               "forward"/"q" - for a forward shift operator
        #               "backward"/"q'" - for a backward shift operator.
        if isa(x, String)
            constructFromString(x)
        else
            if isa(x, Number)
                x = reshape([x],(1,1))
            elseif isa(x, Vector)
                n = length(x)
                x = reshape(x, (n,1))
            end
            n,m = size(x)
            if n < m
                O = Base.zeros(eltype(x), m-n, m)
                x = [x; O]
                n = m
            elseif m < n
                O = Base.zeros(eltype(x), n, n-m)
                x = [x O]
                m = n
            end
            if n > 1
                if iszero(x[end,:]) && iszero(x[:,end])
                    x = x[1:end-1, 1:end-1]
                    ShiftOperator(x)
                else
                    new(x)
                end
            else
                new(x)
            end
        end
    end
end

function constructFromString(str :: String)
    # CONSTRUCTFROMSTRING used by the ShiftOperator constructor if the input
    # is a string
    if str == "q" || str == "forward"
        ShiftOperator([0 1])
    elseif str == "q'" || str == "backward"
        ShiftOperator([0 1]')
    else
        error("Invalid string")
    end
end

function getCoefficients(op::ShiftOperator)
    # GETCOEFFICIENTS Returns the shiftCoefficients of a shift operator.
    # GETCOEFFICIENTS(op) returns the matrix M of the operator op
    # where op = [1 q q²...]M[1 q' (q')^2 ...]'.
    return op.shiftMatrix
end

function op2str(op::ShiftOperator)
    # OP2STR Returns the string representation of a shift operator.
    operatorString = ""

    # Check if it is the zero operator
    z = ShiftOperator(0)
    if(op == z)
        operatorString = "0"
        return operatorString
    end


    if (op.shiftMatrix isa AbstractVector)
        n = length(op.shiftMatrix)
        m = 1
    else
        n,m = size(op.shiftMatrix)
    end

    for i = 1:n
        for j = 1:m
            coefficient = op.shiftMatrix[i,j]

            term = ""
            # If coefficient is non-zero, build the corresponding term
            if (coefficient != 0)
                powerQ = j - 1
                powerQstar = i - 1

                # Add coefficient to term if it should be included
                if (coefficient != 1)
                    term = coefficient;
                elseif (i == 1 && j == 1)
                    term = coefficient
                end

                # Add backward shift operators
                if (powerQstar == 1)
                    term = string( term, "q'")
                elseif (powerQstar != 0)
                    term = string( term, "(q')^", powerQstar)
                end

                # Add forward shift operators
                if (powerQ == 1)
                    term = string( term, "q")
                elseif (powerQ != 0)
                    term = string( term, "q^", powerQ)
                end

                # Add term term to string
                operatorString = string(operatorString, term,"+")
            end

        end
    end

    operatorString = operatorString[1:end-1]
    return operatorString
end



function Base.show(io::IO, op::ShiftOperator)
# Function of displaying shift operators
    operatorString = op2str(op)
    print(io, operatorString)
end

function Base.:(==)(op1 :: ShiftOperator, op2 :: ShiftOperator)
    # == Checks if two matrices of shift operators are equal
    return op1.shiftMatrix == op2.shiftMatrix
end

function Base.:+( op1 :: ShiftOperator, op2 :: ShiftOperator)
    # Function to add two shift operators

    coefficients1 = op1.shiftMatrix
    coefficients2 = op2.shiftMatrix

    # Make sure dimension of coefficient matrices match
    n1 = size(coefficients1,1) 
    m1 = size(coefficients1,2)
    n2 = size(coefficients2,1)
    m2 = size(coefficients2,2)
    n = max(n1,n2)
    m = max(m1,m2)
    extendedCoeff1 = zeros(n,m)
    extendedCoeff2 = zeros(n,m)
    extendedCoeff1[1:n1,1:m1] .= coefficients1
    extendedCoeff2[1:n2,1:m2] .= coefficients2

    sumOperator = ShiftOperator(extendedCoeff1 + extendedCoeff2)
    return sumOperator
end

function Base.:-( op1 :: ShiftOperator, op2 :: ShiftOperator)
    # Function to add two shift operators

    coefficients1 = op1.shiftMatrix
    coefficients2 = op2.shiftMatrix

    # Make sure dimension of coefficient matrices match
    n1 = size(coefficients1,1) 
    m1 = size(coefficients1,2)
    n2 = size(coefficients2,1)
    m2 = size(coefficients2,2)
    n = max(n1,n2)
    m = max(m1,m2)
    extendedCoeff1 = zeros(n,m)
    extendedCoeff2 = zeros(n,m)
    extendedCoeff1[1:n1,1:m1] = coefficients1
    extendedCoeff2[1:n2,1:m2] = coefficients2

    diffOperator = ShiftOperator(extendedCoeff1 - extendedCoeff2)
    return diffOperator
end

function Base.:-(op :: ShiftOperator)
    # Returns the operator with a shifted sign.
    coefficients = op.shiftMatrix
    negOp = ShiftOperator(-coefficients)
    return negOp
end

function Base.:*( op1 :: ShiftOperator, op2 :: ShiftOperator)
    # Function for multiplyig the operators op1 and op2

    # Check if an operator is zero
    z = ShiftOperator(0)
    if (op1 == z || op2 == z)
        return ShiftOperator(0)
    end

    coefficients1 = op1.shiftMatrix
    coefficients2 = op2.shiftMatrix

    n1 = size(coefficients1,1)
    m1 = size(coefficients1,2)
    n2 = size(coefficients2,1)
    m2 = size(coefficients2,2)

    # Initialising the product to be 0
    productCoefficients = 0;

    # Going through all terms in the two operators and multiplying them
    for i = 1:n1
        for j = 1:m1
            coeff1 = coefficients1[i,j]
            if (coeff1 != 0)    #Only continue if the term from the 1st operator isn't 0
                # Multiplyng the term of the 1st operator with all terms of the 2nd
                for k = 1:n2
                    for l = 1:m2
                        coeff2 = coefficients2[k,l]
                        termProductCoefficient = coeff1*coeff2

                        # The multiplication between the two terms on the form
                        # (q^*)^(i-1) q^(j-1) * (q^*)^(k-1) q^(l-1)
                        if (j >= k )
                            # This means that all the (q^*)^(k-1) are cancelled
                            # by the q^(j-1) term. The rest of the q^(j-1) term
                            # is added to the q^(l-1) term
                            
                            # Resize approriatly
                            nOld = size(productCoefficients,1)
                            mOld = size(productCoefficients,2)
                            n = max(nOld, i)
                            m = max(mOld, j+l-k)
                            productCoefficientsNew = zeros(n,m)
                            productCoefficientsNew[1:nOld, 1:mOld] .= productCoefficients
                            productCoefficients = productCoefficientsNew

                            # Add the coefficient product in the spot corresponding
                            # to (q^*)^(i-1) q^(j+l-k)
                            productCoefficients[i, j + l - k] += termProductCoefficient

                        else
                            # This means that all q^(j-1) are cancelled by the
                            # (q^'*)^(k-1) term. The rest of the (q^*)^(k-1) term
                            # is added to the (q^*)^(i-1) term

                            # Resize approriatly
                            nOld = size(productCoefficients,1)
                            mOld = size(productCoefficients,2)
                            n = max(nOld, i+k-j)
                            m = max(mOld, l)
                            productCoefficientsNew = zeros(n,m)
                            productCoefficientsNew[1:nOld, 1:mOld] .= productCoefficients
                            productCoefficients = productCoefficientsNew

                            # Add the coefficient product in the spot corresponding
                            # to (q^*)^(i-1) q^(j+l-k)
                            productCoefficients[i+k-j, l] += termProductCoefficient
                        end
                    end
                end
            end
        end
    end
    return ShiftOperator(productCoefficients)
end

function Base.:*(opMatrix1 :: Matrix{ShiftOperator}, opMatrix2 :: Matrix{ShiftOperator})
    # Multiplies two matrices of shift operators
    n, check1 = size(opMatrix1)
    check2, m = size(opMatrix2)

    if( check1 != check2)
        throw(DimensionMismatch())
    end

    product = Matrix{ShiftOperator}(undef, n, m)
    for i = 1:n
        for j = 1:m
            product[i,j] = ShiftOperator(0)
            for k = 1:check1
                product[i,j] += opMatrix1[i,k]*opMatrix2[k,j]
            end
        end
    end

    return product
end

function Base.:*(matrix :: Matrix, opMatrix :: Matrix{ShiftOperator})
    # Multiplies two matrices, where the right is a matrix of shift operators
    n, check1 = size(matrix)
    check2, m = size(opMatrix)

    if( check1 != check2)
        throw(DimensionMismatch())
    end

    product = Matrix{ShiftOperator}(undef, n, m)
    for i = 1:n
        for j = 1:m
            product[i,j] = ShiftOperator(0)
            for k = 1:check1
                product[i,j] += matrix1[i,k]*opMatrix[k,j]
            end
        end
    end

    return product
end

function Base.:*(opMatrix :: Matrix{ShiftOperator}, matrix :: Matrix)
    # Multiplies two matrices, where the left is a matrix of shift operators
    n, check1 = size(opMatrix)
    check2, m = size(matrix)

    if( check1 != check2)
        throw(DimensionMismatch())
    end

    product = Matrix{ShiftOperator}(undef, n, m)
    for i = 1:n
        for j = 1:m
            product[i,j] = ShiftOperator(0)
            for k = 1:check1
                product[i,j] += opMatrix1[i,k]*matrix[k,j]
            end
        end
    end

    return product
end

function Base.:*(opMatrix :: Matrix{ShiftOperator}, opVector :: Vector{ShiftOperator})
    # Multiplies a matrix of shift operatrs with a vector of shift operators
    opVecMatrix = reshape(opVector, length(opVector), 1)
    return opMatrix*opVecMatrix
end

function Base.:*(opVector :: Vector{ShiftOperator}, opMatrix :: Matrix{ShiftOperator})
    # Multiplies a vector of shift operatrs with a matrix of shift operators
    opVecMatrix = reshape(opVector, length(opVector), 1)
    return opVecMatrix*opMatrix
end


function Base.:*(opVector :: Vector, op :: ShiftOperator)
    # Multiplies a vector of shift operatrs with a shift operator
    n = length(opVector)

    product = Vector{ShiftOperator}(undef, n)
    for i = 1:n
        product[i] = opVector[i]*op
    end

    return product
end

function Base.:*(op :: ShiftOperator, opVector :: Vector)
    # Multiplies a shift operator with a vector of shift operators
    n = length(opVector)

    product = Vector{ShiftOperator}(undef, n)
    for i = 1:n
        product[i] = op*opVector[i]
    end

    return product
end
    
function Base.:*( opMatrix :: Matrix, op :: ShiftOperator)
    # Multiplies a matrix of shift operatrs with a shift operator
    n,m = size(opMatrix)

    product = Matrix{ShiftOperator}(undef, n, m)
    for i = 1:n
        for j= 1:m
            product[i,j] = opMatrix[i,j]*op
        end
    end

    return product
end

function Base.:*( a :: Number, op :: ShiftOperator)
    # Multiplies a scalar with a shift operator
    aOp = ShiftOperator([a])
    return aOp*op
end

function Base.:*( op :: ShiftOperator, a :: Number)
    # Multiplies a shift operator with a scaler
    aOp = ShiftOperator([a])
    return op*aOp
end

function Base.:*(op :: ShiftOperator, opMatrix :: Matrix{ShiftOperator})
    # Multiplies a shift operator with a matrix of shift operators
    n,m = size(opMatrix)

    product = Matrix{ShiftOperator}(undef, n, m)
    for i = 1:n
        for j= 1:m
            product[i,j] = op*opMatrix[i,j]
        end
    end

    return product
end

function Base.:*(op :: ShiftOperator, y :: SignalL2)
    # Applies the operator op to a signal y
    yOut = applyOperator(op,y)
    return yOut
end

function Base.:*(opMatrix :: Matrix{ShiftOperator}, yVec :: Vector{SignalL2})
    # Applies a matrix of operators to a vector of signals
    yOutVec = applyOperator(opMatrix, yVec)
    return yOutVec
end

function Base.:^(op :: ShiftOperator, exponent :: Int)
    # Raises a matrix of shift operators to a power.
    resultOperator = ShiftOperator(1)
    for i = 1:exponent
        resultOperator = resultOperator*op
    end
    return resultOperator
end

function inv(op::ShiftOperator)
    # Inverse for the R-infinty class
    
    coefficients = op.shiftMatrix
    n = size(coefficients,1)
    m = size(coefficients,2)

    # Making sure there are no terms (q^*)^k q^l with k != l
    # If there are, I do not know how to define the inverse
    for i = 1:n
        for j = 1:m
            if (i != j)
                shouldBeZero = coefficients[i,j]
                if (shouldBeZero != 0)
                    throw(DomainError(op, "operator must be in R-infinity class"))
                end
            end
        end
    end

    # Pick out diagonal elements
    diagonalElements = Vector{Float64}(undef, min(n,m))
    for i = 1:n
        for j = 1:m
            if (i == j)
                diagonalElements[i] = coefficients[i,i]
            end
        end
    end

    # Creat partial sums of diagonal elemetns
    N = length(diagonalElements)
    partialSum = Vector{Float64}(undef, N)
    partialSum[1] = diagonalElements[1]
    for i = 2:N
        partialSum[i] = partialSum[i-1] + diagonalElements[i]
    end

    # Check if the operator is invertable
    notInvertable = any( x -> x == 0, partialSum)
    if (notInvertable)
        throw(DomainError(op, "operator is not invertable"))
    end

    # Computing the (q^*)^k q^k coefficients for the inverse
    inverseDiagonalCoefficients = Vector{Float64}(undef,N)
    inverseDiagonalCoefficients[1] = 1/partialSum[1]
    for i = 2:N
        inverseDiagonalCoefficients[i] = 1/partialSum[i] - 1/partialSum[i-1]
    end

    # Build Coefficient matrix
    inverseCoefficients = zeros(Float64, N, N)
    for i = 1:N
        inverseCoefficients[i,i] = inverseDiagonalCoefficients[i]
    end

    return ShiftOperator(inverseCoefficients)

end

function Base.sqrt( op :: ShiftOperator )
    # Square root for the R-infinty class
    
    coefficients = op.shiftMatrix
    n = size(coefficients,1)
    m = size(coefficients,2)

    # Making sure there are no terms (q^*)^k q^l with k != l
    # If there are, I do not know how to define a square root
    for i = 1:n
        for j = 1:m
            if (i != j)
                shouldBeZero = coefficients[i,j]
                if (shouldBeZero != 0)
                    throw(DomainError(op, "operator must be in R-infinity class"))
                end
            end
        end
    end


    # Pick out diagonal elements
    diagonalElements = Vector{Float64}(undef, min(n,m))
    for i = 1:n
        for j = 1:m
            if (i == j)
                diagonalElements[i] = coefficients[i,i]
            end
        end
    end

    # Creat partial sums of diagonal elemetns
    N = length(diagonalElements)
    partialSum = Vector{Float64}(undef, N)
    partialSum[1] = diagonalElements[1]
    for i = 2:N
        partialSum[i] = partialSum[i-1] + diagonalElements[i]
    end


    # Check if the operator is positive definite
    notPosDefinite = any( x -> x < 0, partialSum)
    if (notPosDefinite)
        throw(DomainError(op, "operator must be positive semi-definite"))
    end


    # Computing the (q^*)^k q^k coefficients for the square root
    squareRootDiagonalCoefficients = Vector{Float64}(undef,N)
    squareRootDiagonalCoefficients[1] = sqrt(partialSum[1])
    for i = 2:N
        squareRootDiagonalCoefficients[i] = sqrt(partialSum[i]) - sqrt(partialSum[i-1])
    end

    # Build Coefficient matrix
    squareRootCoefficients = zeros(Float64, N, N)
    for i = 1:N
        squareRootCoefficients[i,i] = squareRootDiagonalCoefficients[i]
    end

    return ShiftOperator(squareRootCoefficients)
end

function Base.adjoint(op :: ShiftOperator)
    # Finds the adjoint of a shift operator
    coefficients = op.shiftMatrix
    transposedOp = ShiftOperator(coefficients')
    return transposedOp
end

function round(op :: ShiftOperator)
    # Rounds the coefficient in a shift operator
    coeff = op.shiftMatrix
    roundedCoeff = Base.round.(coeff)
    return ShiftOperator(roundedCoeff)
end

#------------------- Matrix operatorion --------------------------------


function shiftopZeros(::Type{ShiftOperator}, n ::Int, m :: Int)
    # Gives a shift operator matrix of size n times m with the zero operator
    z = ShiftOperator(0)
    M = Matrix{ShiftOperator}(undef, n, m)
    for i = 1:n
        for j = 1:m
            M[i,j] = z
        end
    end
    return M
end

function shiftopZeros(::Type{ShiftOperator}, n :: Int)
    # Gives a shift operator matrix of size n times n with the zero operator
    return shiftopZeros(ShiftOperator, n, n)
end

function shiftopI(::Type{ShiftOperator}, n :: Int)
    # Gives a shift operator identity matrix of size n
    id = ShiftOperator(1)
    identity = shiftopZeros(ShiftOperator, n)
    for i = 1:n
        identity[i,i] = id
    end
    return identity
end

function maxOrder( M :: Matrix{ShiftOperator})
    # Finds the maximum order of the operators in a matrix of shift operators
    # The maximal order of a shift operator is the maximal power of
    # the shift operators making up the operator.
    n,m = size(M)
    d = 0
    for i = 1:n
        for j= 1:m
            op = M[i,j]
            coeff = getCoefficients(op)
            opOrder = size(coeff, 1) - 1

            d = max(opOrder, d)
        end
    end
    return d
end

function Base.adjoint( opMatrix :: Matrix{ShiftOperator})
    # Finds the adjoint to a matrix of shift operators
    n,m = size(opMatrix)

    transposedOpMatrix = Matrix{ShiftOperator}(undef,m,n)
    for i = 1:n
        for j = 1:m
            op = opMatrix[i,j]
            transposedOpMatrix[j,i] = op' 
        end
    end
    return transposedOpMatrix
end

function chol(M :: Matrix{ShiftOperator})
    # CHOL Performs a Cholesky-likefactorisation of matrices corrseponding to tree graph.
    # CHOL(M) returns the lower triangular matrix of shift operators
    # L such that
    #       M'*M = L*L'
    # M must be a matrix of shift operators with the structure of an
    # incidence matrix of a tree graph. The order of the column in M
    # must correspond to a perfect elimination ordering, i.e. the
    # columns must be order so that each iteration removes a leaf
    # (an edge that is only connected through one node)

    n,m = size(M)

    if (n-1 != m)
        throw(DomainError(M, "there must be exactly one more row the column.\notherwise the matrix cannot represent a tree"))
    end

    # Base case
    if m == 1
        N = M'*M
        N = N[1]
        L = sqrt(N)
        return L
    end

    # Look if vertex permutation is needed
    z = ShiftOperator(0)
    id = ShiftOperator(1)
    if M[1,1] == z
        candidateConnections = M[2:end, 1]
        connection = findFirstElement(candidateConnections) + 1
        P = shiftopI(ShiftOperator, n)
        P[1,1] = z
        P[1,connection] = id
        P[connection,1] = id
        P[connection,connection] = z
        M = P*M
    end
    if M[2,1] == z
        candidateConnections = M[3:end,1]
        connection = findFirstElement(candidateConnections) +2
        P = shiftopI(ShiftOperator, n)
        P[2,2] = z
        P[2, connection] = id
        P[connection, 2] = id
        P[connection, connection] = z
        M = P*M
    end
    M12 = M[1, 2:end]
    if !onlyZeros(M12)
        P = shiftopI(ShiftOperator, n)
        P[1,1] = z
        P[1,2] = id
        P[2,1] = id
        P[2,2] = z
        M = P*M
    end

    # Extract blocks
    M11 = M[1,1]
    M12 = M[1, 2:end]
    M21 = M[2,1]
    M22 = M[2:2, 2:end]
    M31 = M[3:end, 1:1]
    M32 = M[3:end, 2:end]

    # Ensure that M has the correct structure
    if !onlyZeros(M12) || !onlyZeros(M31)
        throw(DomainError(M, "M must have the structure of a tree with an elimination ordering of the edges"))
    end

    # Compute the necessary blocks in the product N = M'*M
    N11 = M11'*M11 + M21'*M21
    N21 = M22'*M21

    # Compute reduced M
    Mred =[sqrt(id - M21*inv(N11)*M21')*M22; M32]

    # Induction step
    Lred = chol(Mred)

    # Get the factor
    O = shiftopZeros(ShiftOperator, 1, m-1)
    N11Sqrt = sqrt(N11)
    L = [N11Sqrt O; N21*inv(N11Sqrt) Lred]

    return L
end

function monomialCoefficients(M :: Matrix{ShiftOperator}, I, J)
    # MONOMIALCOEFFICIENTS Picks out the coefficients for one type of monomial.
    # MONOMIALCOEFFICIENTS(M,i,j) picks the real valued matrix of
    # coefficients of the  monomial (q')^(j-1)q^(i-1) from M.
    #
    # Example:
    # q = ShiftOperator("q")
    # M = [2*q q'; q+q' 3*q^2]
    # Mq = monomialCoefficients(M,1,2)   Picks out the matrix [2 0;1 0]

    n,m = size(M)
    Mpart = zeros(n, m)
    for i = 1:n
        for j = 1:m
            op = M[i,j]
            coeff = op.shiftMatrix
            nOp, mOp = size(coeff)
            if I > nOp || J > mOp
                Mpart[i,j] = 0
            else
                Mpart[i,j] = coeff[I,J]
            end
        end
    end
    return Mpart
end

function findFirstElement(V :: Vector{ShiftOperator})
    # Finds the first non-zeros operator in a vector of shift operators
    n = length(V)
    z = ShiftOperator(0)
    for i = 1:n
        vi = V[i]
        if vi != z
            return i
        end
    end
    return []
end

function onlyZeros(M :: Matrix{ShiftOperator})
    # Checks if the matrix of shift operators only contains the zero operator
    n,m = size(M)
    O = shiftopZeros(ShiftOperator, n, m)
    return O == M
end

function onlyZeros(V = Vector{ShiftOperator})
    # Checks if the vector of shift operators only contains the zero operator
    n = length(V)
    Vm = reshape(V, n, 1)
    O = shiftopZeros(ShiftOperator, n, 1)
    return O == Vm
end

function applyOperator(op::ShiftOperator, y::SignalL2)
    # Applies an operator to a signal
    if op == ShiftOperator(0)
        newY = SignalL2(reshape([0],(1,1)),reshape([0],(1,1)),[0])
        return newY
    end
    coeffs = getCoefficients(op)
    n,m = size(coeffs)
    noTermAdded = true
    for i = 1:n
        for j = 1:m
            cij = coeffs[i,j]
            if (cij != 0)
                yij = applyMonomialOperator(cij, i, j, y)
                if noTermAdded
                    newY = yij
                    noTermAdded = false
                else
                    newY = newY + yij;
                end
            end
        end
    end
    return newY
end

function applyOperator(M :: Matrix{ShiftOperator}, yVec :: Vector{SignalL2})
    # Applies a matrix of operators to a vector of signals
    n,m = size(M)
    nVec, = size(yVec)
    if nVec != m
        error("Number of columns in M must be the same as the number of rows in yVec")
    end
    yNew = [];
    for i = 1:n
        yi = SignalL2(reshape([0],(1,1)),reshape([0],(1,1)),[0])
        for j = 1:m
            yi = yi + applyOperator(M[i,j], yVec[j])
        end
        yNew = [yNew; yi];
    end
    return yNew
end
