using LinearAlgebra
struct SignalL2
    C       # Output vector
    A       # State matrix
    x0      # Initial condition

    # SIGNALL2 
    # This class represents a set of signals in l2+ that can be represented
    # by a tripple such that the signal
    #   w = (C x_0, C A x_0, C A^2 x_0, ...)
    # is represented by the tripple (C, A, x0). The appropriate sizes are
    #   C  - 1 x n (row vector)
    #   A  - n x n (square matrix)
    #   x0 - n x 1 (column vector)
    #
    # The class can be used together with the shiftOperator class. This
    # class can be used when applying such an operator to a signal.
    #
    # Examples:
    # q = shiftOperator('q');
    # op = q^2 + q';
    # syms C A x0;
    # y = signalL2(C, A, x0);
    # yOut = applyOperator(op, y);  % Alterative syntax: op*y

    function SignalL2(C::Matrix, A::Matrix, x0::Vector)
        # Constructs an instance of this class

        nC, mC = size(C);
        
        # Check size A
        nA, mA = size(A)
        if nA != mA
            error("A must be a square matrix")
        elseif mA != mC
            error("A and C must have the same number of columns")
        end
        
        nx, = size(x0)
        if nx != nA
            error("x0 must have the same number of rows as A")
        end
        C, A, x0 = dimensionReduction(C, A, x0)
        new(C, A, x0)
    end
end

function getC(y::SignalL2)
    # GETC returns the output matrix
    # It gives the vector C of the signal
    #   w = (C x0, C A x0, C A^2 x0, ...)
    return y.C
end

function getA(y::SignalL2)
    # GETA returns the state matrix
    # It gives the matrix A of the signal
    #   w = (C x0, C A x0, C A^2 x0, ...)
    return y.A
end

function getX0(y::SignalL2)
    # GETX0 returns the initial condition
    # It gives the vector x0 of the signal
    #   w = (C x0, C A x0, C A^2 x0, ...)
    return y.x0;
end

function firstSamples(y::SignalL2, k::Int)
    # Gives a strin showing the first k entries of a signal
    signalString ="("
    for i = 0:k-1
        sample = reshape(y.C * y.A^i * y.x0, (1))
        signalString = string(signalString, sample, ", ")
    end
    signalString = string(signalString, "...)")
    return signalString
end

function Base.show(io::IO, y::SignalL2)
    # Displays signals
    k = 4
    signalString = firstSamples(y, k)
    print(io, signalString)
end

function Base.:+(y1::SignalL2, y2::SignalL2)
    # Addition of signals
    Csum = [y1.C y2.C]

    n1,m1 = size(y1.A)
    n2,m2 = size(y2.A)
    Asum = [y1.A zeros(n1,n2); zeros(n2,n1) y2.A]

    x0sum = [y1.x0; y2.x0]

    signalSum = SignalL2(Csum, Asum, x0sum)
    return signalSum
end

function Base.:-(y1::SignalL2, y2::SignalL2)
    # Subtraction of signals
    y2Neg = SignalL2(-y2.C, y2.A, y2.x0)
    signalDiff = y1 + y2Neg
    return signalDiff
end

function applyMonomialOperator(c::Number, i::Int, j::Int, y::SignalL2)
    # Applies an operator c*(q')^i*q^j to a signal y
    if i < j
        newC = c*y.C*y.A^(j-i)
        newA = y.A
        newX0 = y.x0
    else
        n,m = size(y.A)
        newC = [zeros(1,(i-1)*n) c*y.C*y.A^(j-1)]
        newA = zeros(n*i,n*i)
        newA[1:n,1:n] = y.A
        if i > 1
            Idn = Matrix{Int}(I, n*(i-1), n*(i-1))
            newA[n+1:end,1:end-n] = Idn
        end
        newX0 = zeros(n*i)
        newX0[1:n] = y.x0
    end
    newY = SignalL2(newC, newA, newX0)
    return newY
end

function dimensionReduction(C::Matrix, A::Matrix, x0::Vector)
    # Reduces the dimension for the minimum necessary dimension for observability
    n,m = size(A)
    Wo = C
    for i = 2:n
        Wo = [Wo; C*A^(i-1)]
    end
    n = rank(Wo)
    if n == 0
        Cred = reshape([0], (1,1))
        Ared = reshape([0], (1,1))
        x0red = [0]
        return Cred, Ared, x0red
    end
    
    u,s,v = svd(Wo)
    Ared = v'*A*v
    Ared = Ared[1:n, 1:n]
    Cred = C*v
    Cred = Cred[1, 1:n]
    Cred = reshape(Cred,(1,n))
    x0red = v'*x0
    x0red = x0red[1:n, 1]
    return Cred, Ared, x0red
end
