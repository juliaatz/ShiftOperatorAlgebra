function [K1,K2] = factor_lqr_formula(A, C, M, r)
    L = chol(M);
    L0 = monomialCoefficients(L,1,1);
    L1 = monomialCoefficients(L,1,2);
    
    Mstar0 = monomialCoefficients(M',1,1);
    Mstar1 = monomialCoefficients(M',1,2);
    
    q = shiftOperator('q');
    if (L0 + L1*q ~= L) || (Mstar0 + Mstar1*q ~= M')
        error("The formula does not hold. Try factor_lqr instead")
    end
    
    K1 = (L0 + L1*r)*L0';
    K2 = (Mstar0 + Mstar1*r)*C*r*A;
end