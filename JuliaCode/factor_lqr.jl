function factor_lqr(A :: Matrix, C :: Matrix, M :: Matrix{ShiftOperator}, r :: Number)
    M0 = monomialCoefficients(M,1,1)
    M1 = monomialCoefficients(M,2,1)

    q = ShiftOperator("q")
    id = ShiftOperator(1)

    if M0*id + M1*q' != M
        throw(DomainError(M, "M must be expressable as M0 + M1*q"))
    end

    n,check = size(A)
    if n != check
        throw(DomainError(A,"A must be squar"))
    end
    p,m = size(M)

    # Find matrix representation of L
    L = 0
    try
        L = chol(M)
    catch
        throw(DomainError(M, "M could not be factorised with chol"))
    end
    d = maxOrder(M)
    Lm = matrixRepresentation(L, d)

    # Initial system, corresponding to q y_init
    A1 = r*A
    C1 = r*C*A

    # Second system, corresponding to M'
    Ip = zeros(p,p)
    for i = 1:p
        Ip[i,i] = 1
    end
    B2 = [zeros((d-1)*p,p); Ip]

    A2 = zeros(d*p,d*p)
    C2 = zeros(m, d*p)
    for i = 0:(d-1)
        if i < (d-1)
            A2[p*i+1:p*(i+1), (i+1)*p+1:(i+2)*p] = Ip
        end

        Mi = monomialCoefficients(M', 1, i+1)
        C2[:, p*i+1:(i+1)*p] = Mi
    end
    Md = monomialCoefficients(M', 1, d+1)
    D2 = Md

    # Feed systems into one another
    Ac = [A1 zeros(n,d*p); B2*C1 A2]
    Cc = [D2*C1 C2]

    l,l = size(Lm)
    N = Int8(l/m - 1)
    wsExtendedState = zeros((N+1)*m, n+p*d)
    for i = 0:N
        wsExtendedState[m*i+1:(i+1)*m, :] = Cc*Ac^(i+d)
    end
    ws = wsExtendedState[:, 1:n]

    # Build factors
    w0 = ws[1:m, :]
    wFactor = zeros(l,m)
    for i = 0:N
        wi = ws[m*i+1:(i+1)*m,:]

        factor = wi/w0
        # Check accuracy of factor
        diff = wi - factor*w0
        errorF = sum(abs.(diff))
        if errorF > 0.001*length(wi)
            error("unexpected error when factoring")
        end

        wFactor[m*i+1:(i+1)*m, :] = factor
    end
    
    # Build factors
    K2 = w0
    nuF = Lm\wFactor
    vF = Lm'\nuF
    K1Inv = vF[1:m,:]
    K1 = Base.inv(K1Inv)
    return K1, K2
end

function matrixRepresentation(L, maxDelay)
    n,m = size(L)
    if n != m
        throw(DomainError(L, "L must be square"))
    end

    # Initialise large matrix
    Lm = Base.zeros(maxDelay*n^2 + maxDelay*n, maxDelay*n^2 + maxDelay*n)
    for i = 0:maxDelay
        Li = monomialCoefficients(L,1,1+i)
        LiMinus = monomialCoefficients(L,i+1,1)
        for j = 1:(n+1)*maxDelay-i
            # Adding upper triangular part
            startRowU = 1 + n*(j-1)
            startColumnU = 1 + n*i + n*(j-1)
            Lm[startRowU:startRowU + (n-1), startColumnU: startColumnU + (n-1)] = Li
            Li = Li + monomialCoefficients(L, 1+j, i+j+1)

            # Adding lower triangular part
            if i != 0
                startRowL = startColumnU
                startColumnL = startRowU
                Lm[startRowL:startRowL+(n-1), startColumnL: startColumnL + (n-1)] = LiMinus
                LiMinus = LiMinus + monomialCoefficients(L,i+j+1,1+j)
            end
        end
    end
    return Lm
end
