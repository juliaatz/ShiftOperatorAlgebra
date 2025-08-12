function [K1, K2] = factor_lqr(A, C, M, r)
%FACTOR_LQR Finds the factorised version of the LQR when possible.
%
%[K1, K2] = FACTOR_LQR(A, C, M, r) gives a factorised version of the LQR on
%the form
%       K1*u(0) = K2*x(0),
%when it is possible. For more details of when the factorisation is
%possible SEE HERE.

M0 = monomialCoefficients(M,1,1);
M1 = monomialCoefficients(M,2,1);
q = shiftOperator('q');
Mexpected = M0 + M1*q';
if ~onlyZeros(M - Mexpected)
    errorStruct.message ="Cannot find factorised control law. M must be expressable as M0+M1*q'";
    errorStruct.identifier = "Flqr:M";
    error(errorStruct);
end

n = length(A);
[p, m] = size(M);

% find matrix representation
try
    L = chol(M);
catch
    errorStruct.message = "Cannot find a factorised control law. Matrix M could not be factorised with chol";
    errorStruct.identifier = "Flqr:chol";
    error(errorStruct);
end
d = maxOrder(M);
Lm = matrixRepresentation(L, d);


% Initial system, correcponding to q y_init
A1 = r*A;
C1 = r*C*A;

% Second system, corresponding to M'
B2 = [zeros((d-1)*p,p); eye(p)];
A2 = zeros(d*p);
C2 = zeros(m,d*p);
for i = 0:(d-1)
    if i < (d-1)
        A2(p*i+1:p*(i+1), (i+1)*p+1:(i+2)*p) = eye(p);
    end
    
    Mi = monomialCoefficients(M', 1, i+1);
    C2(:,p*i+1:(i+1)*p) = Mi;
end
Md = monomialCoefficients(M', 1, d+1);
D2 = Md;

% Feed systems into one another
Ac =[A1 zeros(n,d*p); B2*C1 A2];
Cc = [D2*C1 C2];

l = length(Lm);
N = l/m - 1;
wsExtended = zeros((N+1)*m,n+p*d);
for i = 0:N
    wsExtended(m*i+1:(i+1)*m,:) = Cc*Ac^(i+d);
end
ws = wsExtended(:,1:n);

% Try to factor
w0 = ws(1:m,:);
canFactorise = true;
wFactor = zeros(l,m);
for i = 0:N
    wi = ws(m*i+1:(i+1)*m,:);
    
        factor = wi/w0;
        % Check if accurate factor
        diff = wi - factor*w0;
        errorF = sum(abs(diff),'all');
        if errorF > 0.001*numel(wi)
            factor = 0*eye(m);
            canFactorise = false;
        end
    wFactor(m*i+1:(i+1)*m,:) = factor;
end

% Build factors if possible
if canFactorise
    K2 = w0;
    nuF = Lm\wFactor;
    vF = Lm'\nuF;
    K1Inv = vF(1:m,:);
    K1 = inv(K1Inv);
else
    errorStruct.message = "Cannot find factorisation";
    errorStruct.identifier = "Flqr:unexpected";
    error(errorStruct);
end
end

function Lmatrix = matrixRepresentation(L, maxDelay)
% Represents the matrix of operators as a large matrix of real numbers
[n,m] = size(L);
if (n ~= m)
    error("L not square")
end

% Build big matrix
Lmatrix = zeros(maxDelay*n^2 + maxDelay*n);
if (class(getCoefficients(L)) == "sym")
    Lmatrix = sym(Lmatrix);
end
for i = 0:maxDelay
    Li = monomialCoefficients(L,1,i+1);
    LiMinus = monomialCoefficients(L,i+1,1);
    for j = 1:(n+1)*maxDelay-i
        % Adding upper triangular part
        startRowU = 1 + n*(j-1);
        startColumnU = 1 + n*i + n*(j-1);
        Lmatrix(startRowU:startRowU+(n-1), startColumnU: startColumnU + (n-1)) = Li;
        Li = Li + monomialCoefficients(L,1+j,i+j+1);
        
        % Adding lower triangular part
        if i ~= 0
            startRowL = startColumnU;
            startColumnL = startRowU;
            Lmatrix(startRowL:startRowL+(n-1), startColumnL: startColumnL + (n-1)) = LiMinus;
            LiMinus = LiMinus + monomialCoefficients(L,i+j+1,1+j);
        end
        
    end
    
end
if (class(getCoefficients(L)) == "sym")
    Lmatrix = simplify(Lmatrix);
end
end



