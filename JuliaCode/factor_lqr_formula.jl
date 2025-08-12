function factor_lqr_formula(A :: Matrix, C :: Matrix, M :: Matrix{ShiftOperator}, r :: Number)
    L = chol(M)
    L0 = monomialCoefficients(L,1,1)
    L1 = monomialCoefficients(L,1,2)

    Mstar0 = monomialCoefficients(M',1,1)
    Mstar1 = monomialCoefficients(M',1,2)

    q = ShiftOperator("q")
    id = ShiftOperator(1)
    if (L0*id + L1*q != L) || (Mstar0*id + Mstar1*q != M')
        throw(DomainError(M, "the formula does not hold for this M. Try factor_lqr instead"))
    end

    K1 = (L0 + L1*r)*L0'
    K2 = (Mstar0 + Mstar1*r)*C*r*A
    return K1, K2
end

