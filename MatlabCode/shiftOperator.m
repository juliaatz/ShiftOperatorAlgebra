classdef shiftOperator
    % SHIFTOPERATOR
    % This is a class of shift operators and matrices of shift operator on
    % l2+. These operators are built from the forwardshift operator q and
    % its adjoint the backward shift operator q'.
    %
    % Given a signal y in l2+
    %       q : [y1, y2, y3, ...] -> [y2, y3, y4, ...]
    %       q': [y1, y2, y3, ...] -> [0, y1, y2, ...].
    %
    % This class allows you to perform algebraic operations with these
    % obejcts.
    
    properties
        shiftCoefficients              % matrix M for a operate op such that op = [1 q q^2...]M[1 q' (q')^2 ...]'.
    end
    
    methods
        function op = shiftOperator(shiftData)
            %SHIFTOPERATOR Constructs a shift operator from its shift coefficients or a string.
            %SHIFTOPERATOR(shiftData) creates a shift operator.
            %ShiftData can either be a matrix or a string.
            
            %Matrix: shiftData is the matrix M of coefficeients such that
            %the operator is [1 q q^2...]M[1 q' (q')^2 ...]'.
            
            %String: ShiftData is on of the following 4 strings:
            %               "forward"/"q" - for a forward shift operator
            %               "backward"/"q'" - for a backward shift operator.
            shiftsExists = exist('shiftData','var');
            if (~shiftsExists)
                shiftData = 0;
            end
            if isstring(shiftData) || ischar(shiftData)
                if strcmp(shiftData, "forward") || strcmp(shiftData, "q")
                    shiftData = [0 1];
                elseif strcmp(shiftData, "backward") || strcmp(shiftData, "q'")
                    shiftData = [0 1]';
                else
                    errorStruct.message = "Input to constructor must be a matrix of coefficients or one of the following strings: forward, q, backward, q'";
                    errorStruct.identifier ="Constructor:invalidString";
                    error(errorStruct)
                end
            end
            
            % Make quadratic for easier code
            [n,m] = size(shiftData);
            if (n > m)
                shiftData(n,n) = 0;
            elseif (m > n)
                shiftData(m,m) = 0;
            end
            
            % Ensure that the matrix isn't larger than needed
            neededElementFound = false;
            while (~neededElementFound && n > 1)
                lastRow = shiftData(end,:);
                lastColumn = shiftData(:,end);
                
                nbrElementsRow = length(find(lastRow));
                nbrElementsColumn = length(find(lastColumn));
                
                if( nbrElementsColumn + nbrElementsRow == 0)
                    shiftData = shiftData(1:n-1, 1:n-1);
                    n = n-1;
                else
                    neededElementFound = true;
                end
            end
            
                
            op.shiftCoefficients = shiftData;
        end
        
        function disp(op)
            %DISP Displays shift operators and matrices of shift operators.
            %DISP(op) displays op where op can be either a shift operator
            %or a matrix of shift operators.
            [n,m] = size(op);
            if(n == 1 && m == 1)
                operatorRow = op2str(op);
                disp(operatorRow)
            else
                for i = 1:n
                    operatorRow =  "[";
                    for j = 1:m
                        operatorElement = op2str(op(i,j));
                        operatorRow = strcat(operatorRow, operatorElement,',');
                    end
                    operatorRow = convertStringsToChars(operatorRow);
                    operatorRow = operatorRow(1:end-1);
                    operatorRow = strcat(operatorRow, "]");
                    disp(operatorRow)
                end
            end
        end
        
        function coeff = getCoefficients(op)
            %GETCOEFFICIENTS Returns the shiftCoefficients of a shift operator.
            %GETCOEFFICIENTS(op) returns the matrix M of the operator op
            %where op = [1 q q²...]M[1 q' (q')^2 ...]'.
            coeff = op.shiftCoefficients;
        end
            
        function operatorString = op2str(op)
            %OP2STR Returns the string representation of a shift operator.
            operatorChar = '';
            if (isempty(op))
                return
            end
            shiftMatrix = op.shiftCoefficients;
            
            [n,m] = size(shiftMatrix);
            for i = 1:n
                for j = 1:m
                    coefficient = shiftMatrix(i,j);
                    if (abs(coefficient) ~= 0)
                        powerQ = j-1;
                        powerQstar = i-1;
                        
                        monomial = "";
                        if (coefficient ~= 1)
                            if (class(coefficient) == "double")
                                monomial = num2str(coefficient);
                            else
                                monomial = string(coefficient);
                            end
                        elseif (i==1 && j==1)
                            monomial = "1";
                        end
                        
                        if (powerQstar == 1)
                            monomial = strcat(monomial, "q'");
                        elseif (powerQstar ~= 0)
                            monomial = strcat(monomial, "(q')^", int2str(powerQstar));
                        end
                        if (powerQ == 1)
                            monomial = strcat(monomial, 'q');
                        elseif (powerQ ~= 0)
                            monomial = strcat(monomial, 'q^', int2str(powerQ));
                        end
                        operatorChar = strcat(operatorChar,monomial,'+');
                    elseif (sum(abs(shiftMatrix),'all') == 0 )
                        operatorChar = "0";
                        operatorString = convertCharsToStrings(operatorChar);
                        return
                    end
                end
            end
            operatorChar = convertStringsToChars(operatorChar);
            operatorChar = operatorChar(1:end-1);
            operatorString = convertCharsToStrings(operatorChar);
        end
        
        function negOperator = uminus(operators)
            %UMINUS Returns the operator with a shifted sign.
            %UMINUS(op) is the same as -op
            [n,m] = size(operators);
            for i = 1:n
                for j = 1:m
                    op = operators(i,j);
                    negOperator(i,j) = shiftOperator(-op.shiftCoefficients);
                end
            end
        end
        
        function opSum = plus(operator1, operator2)
            %PLUS Sums two matrices of shift operators.
            %PLUS(A,B) is the same as A+B where A and B are matrices of
            %shift operators of the same size.
            [n1, m1] = size(operator1);
            [n2, m2] = size(operator2);
            
            if (n1 ~= n2 || m1 ~= m2)
                errorStruct.message = "Matrices must have the same dimensions";
                errorStruct.identifier = "MatrixAddition:dimensionMismatch";
                error( errorStruct )
            end
            
            for i = 1:n1
                for j = 1:m1
                    op1 = operator1(i,j);
                    op2 = operator2(i,j);
                    opSum(i,j) = plusElement(op1,op2);
                end
            end
        end
       
        function opDiff = minus(operator1, operator2)
            %MINUS Subtraction of matrices of shift operators.
            %MINUS(A,B) is the same as A-B where A and B are matrices of
            %shift operators of the same size.
            [n1, m1] = size(operator1);
            [n2, m2] = size(operator2);
            
            if (n1 ~= n2 || m1 ~= m2)
                errorStruct.message = "Matrices must have the same dimensions";
                errorStruct.identifier = "MatrixSubtraction:dimensionMismatch";
                error( errorStruct )
            end
            
            for i = 1:n1
                for j = 1:m1
                    op1 = operator1(i,j);
                    op2 = operator2(i,j);
                    opDiff(i,j) = minusElement(op1,op2);
                end
            end
        end
        
        function matrixProduct = times(A,B)
            %TIMES Element wise matrix multiplication of shift operators.
            %TIMES(A,B) gives the same result as A.*B where A and B are
            %shift operators of the same size.
            
            %Ensure correct dimenssions
            [n, m] = size(A);
            [nb, mb] = size(B);
            
            if n~=nb || m~=mb
                errorStruct.message = "Matrices must have same dimensions.";
                errorStruct.identifier = "Times:dimensions";
                error(errorStruct)
            end
            
            matrixProduct = shiftOperator.zeros(n,m);
            for i = 1:n
                for j = 1:m
                    matrixProduct(i,j) = etimes(A(i,j), B(i,j));
                end
            end
        end
        
        function matrixProduct = mtimes(A,B)
            %MTIMES Matrix multiplication of matrices of shift operators.
            %MTIMES(A,B) gives the same result as A*B where A and B are
            %matrices of shift operators of approriate size.
            [n, checkA] = size(A);
            [checkB, m] = size(B);
            if(checkA ~= checkB)
                if (n == 1) && (checkA == 1)
                    matrixProduct = operatorMatrixProduct(A,B);
                    return
                elseif (m == 1) && (checkB == 1)
                    matrixProduct = matrixOperatorProduct(A,B);
                    return
                end
                errorStruct.message = "The number of columns in A must be the same as the number of rows in B";
                errorStruct.identifier = "MatrixMultiplication:dimensionMismatch";
                error(errorStruct);
            end
            
            for i = 1:n
                for j = 1:m
                    matrixProduct(i,j) = shiftOperator(0);
                    for k = 1: checkA
                        matrixProduct(i,j) = matrixProduct(i,j) + A(i,k).*B(k,j);
                    end
                end
            end
        end
        
        function transposedOperator = ctranspose(operators)
            %CTRANSPOSE Finds the adjoint of an operator.
            %CTRANSPOSE(op) is the same as op'.
            [n,m] = size(operators);
            
            for i = 1:n
                for j = 1:m
                    op = operators(i,j);
                    newShifts = op.shiftCoefficients;
                    transposedOperator(j,i) = shiftOperator(newShifts');
                end
            end
        end
        
        function inverse = inv(op)
            %INV Finds the inverse of operators in the R-infinity class.
            shiftMatrix = op.shiftCoefficients;
            
   
            if (~isRinfOperator(op))
                errorStruct.message = "The operator you are trying to invert is not in the R-infinity set. The inverse is not known";
                errorStruct.identifier = "Inverse:nonRinf";
                error(errorStruct);
            end
            
            diagPart = diag(shiftMatrix);
            n = length(diagPart);
            s = zeros(1,n);
            if (class(diagPart) == "sym")
                s = sym(s);
            end
            s(1) = diagPart(1);
            for i = 2:n
                s(i) = s(i-1) + diagPart(i);
            end
            
            % Check if invertable
            zeroSums = find(s==0);
            if(~isempty(zeroSums))
                errorStruct.message = "The operator is not invertable";
                errorStruct.identifier = "Inverse:nonInvertable";
                error(errorStruct)
            end
            
            inverseDiag = zeros(1,n);
            if (class(diagPart) == "sym")
                inverseDiag = sym(inverseDiag);
            end
            inverseDiag(1) = 1/s(1);
            for i = 2:n
                inverseDiag(i) = 1/s(i) - 1/s(i-1);
            end
            if (class(inverseDiag) == "sym")
                inverseDiag = simplify(inverseDiag);
            end
            inverse = shiftOperator(diag(inverseDiag));
        end
        
        function powerOperator = power(operators,a)
            %POWER Element wise power of matrix of shift operators.
            %POWER(M,a) is the same as M.^a where M is a matrix of shift
            %operators and a is a natural number.
            if (class(a) == "shiftOperator")
                errorStruct.message = "Operation is not defined. The exponent has to be a natural number";
                errorStruct.identifier ="Power:operator";
                error(errorStruct)
            end
            if (class(a) =="sym")
                errorStruct.message = "Operation is not defined. The exponent has to be a natural number";
                errorStruct.identifier = "Power:symbolic";
                error(errorStruct)
            end
            if ( a < 0 )
                errorStruct.message = "Operation is not defined. The exponent has to be a natural number";
                errorStruct.identifier ="Power:negative";
                error(errorStruct)
            end
            [n,m] = size(a);
            if n > 1 || m > 1
                errorStruct.message = "Operation is not defined. The exponent has to be a natural number";
                errorStruct.identifier = "Power:matrix";
                error(errorStruct)
            end
            
            [n,m] = size(operators);
            
            for i = 1:n
                for j= 1:m
                    op = operators(i,j);
                    
                    product = shiftOperator(1);
                    for k = 1:a
                        product = product*op;
                    end
                    powerOperator(i,j) = product;
                end
            end
        end
        
        function product = mpower(M,a)
            %MPOWER Raises a matrix of shift operators to a power.
            %MPOWER(M,a) is the same as M^a where M is a matrix of shift
            %operators and a is a natural number.
            if (class(a) == "shiftOperator")
                errorStruct.message = "Operation is not defined. The exponent has to be a natural number";
                errorStruct.identifier ="Mpower:operator";
                error(errorStruct)
            end
            if (class(a) =="sym")
                errorStruct.message = "Operation is not defined. The exponent has to be a natural number";
                errorStruct.identifier = "Mpower:symbolic";
                error(errorStruct)
            end
            if ( a < 0 )
                errorStruct.message = "Operation is not defined. The exponent has to be a natural number";
                errorStruct.identifier ="Mpower:negative";
                error(errorStruct)
            end
            [n,m] = size(a);
            if n > 1 || m > 1
                errorStruct.message = "Operation is not defined. The exponent has to be a natural number";
                errorStruct.identifier = "Mpower:matrix";
                error(errorStruct)
            end
            
            [n,m] = size(M);
            if n~=m
                errorStruct.message = "Matrix must be square";
                errorStruct.identifier = "Mpower:dimensions";
                error(errorStruct)
            end
            
            product = shiftOperator.eye(n);
            for i = 1:a
                product = product*M;
            end
        end
        
        function root = sqrt(op)
            %SQRT Computes the square roots of shift operators in the R-infinity class.
            if(~isRinfOperator(op))
                errorStruct.message = "The operator must be in the R-infinity set";
                errorStruct.identifier = "Sqrt:nonRinf";
                error(errorStruct)
            end
            coefficients = op.shiftCoefficients;
            
            ds = diag(coefficients);
            s = cumsum(ds);
            
            if ( sum( double(sign(s)) < 0 ) > 0 )
                errorStruct.message = "The partial sums are negative. No square root is defined";
                errorStruct.identifier = "Sqrt:negative";
                error(errorStruct)
            end
            
            t = sqrt(s);
            
            n = length(t);
            rootCoef = zeros(1,n);
            if (class(t) == "sym")
                rootCoef = sym(rootCoef);
            end
            
            rootCoef(1) = t(1);
            
            for i = 2:n
                rootCoef(i) = t(i) - t(i-1);
            end
            
            if (class(t) == "sym")
                rootCoef = simplify(rootCoef);
            end
            
            
            root = shiftOperator(diag(rootCoef));
        end
        
        function bool = isRinfOperator(op)
            %ISRINFOPERATOR Checks if a shift operator is in the R-infinity class.
            %ISRINFOPERATOR(op) returns true if op is in R-infinity
            shiftMatrix = op.shiftCoefficients;
            
            offDiagonal = shiftMatrix - diag(diag(shiftMatrix));
            if (sum(abs(offDiagonal)) == 0)
                bool = true;
            else
                bool = false;
            end
            
        end
        
        function L = chol(M)
            %CHOL Performs a Cholesky-likefactorisation of matrices corrseponding to tree graph.
            %CHOL(M) returns the lower triangular matrix of shift operators
            %L such that
            %       M'*M = L*L'
            %M must be a matrix of shift operators with the structure of an
            %incidence matrix of a tree graph. The order of the column in M
            %must correspond to a perfect elimination ordering, i.e. the
            %columns must be order so that each iteration removes a leaf
            %(an edge that is only connected through one node)
 
            
            [n,m] = size(M);
            
            % Ensure that M could correspond to a tree graph
            if ~(n-1 == m)
                errorStruct.messager = "M must have exactly one more row than column. Otherwise can not represent a tree";
                errorStruct.identifier = "Chol:notTree";
                error(errorStruct);
            end
            
            % Base case
            if m == 1
                L = sqrt(M'*M);
                return
            end
            
            % Look if vertex permutation is needed
            if isZeroOperator(M(1,1))
                candidateConnections = M(2:end,1);
                connectionVec = find(candidateConnections) + 1;
                connection = connectionVec(1);
                P = shiftOperator.eye(n);
                P(1,1) = shiftOperator(0);
                P(1, connection) = shiftOperator(1);
                P(connection, 1) = shiftOperator(1);
                P(connection, connection) = shiftOperator(0);
                M = P*M;
            end
            if isZeroOperator(M(2,1))
                candidateConnections = M(3:end,1);
                connectionVec = find(candidateConnections) + 2;
                connection = connectionVec(1);
                P = shiftOperator.eye(n);
                P(2,2) = shiftOperator(0);
                P(2, connection) = shiftOperator(1);
                P(connection, 2) = shiftOperator(1);
                P(connection, connection) = shiftOperator(0);
                M = P*M;
            end
            if ~onlyZeros(M(1, 2:end)) % Switches order of leef nodes if necessary
                P = shiftOperator.eye(n);
                P(1,1) = shiftOperator(0);
                P(1, 2) = shiftOperator(1);
                P(2, 1) = shiftOperator(1);
                P(2, 2) = shiftOperator(0);
                M = P*M;
            end
            
            
            
            % Extract the blocks of M
            M11 = M(1,1);
            M12 = M(1, 2:end);
            M21 = M(2,1);
            M22 = M(2, 2:end);
            M31 = M(3:end, 1);
            M32 = M(3:end,2:end);
            
            % Ensure that M has the correct structure
            if (~onlyZeros(M12)) || (~onlyZeros(M31))
                errorStruct.message = "The matrix M must have the structure of a tree, with an elimination ordering of the edges";
                errorStruct.identifier = "Chol:notLeaf";
                error(errorStruct);
            end
            
            % Compute necessary blocks in product N = M'*M
            N11 = M11'*M11 + M21'*M21;
            N21 = M22'*M21;
            
            % Compute reduced M
            one = shiftOperator(1);
            Mred = [sqrt(one - M21*inv(N11)*M21')*M22; M32];
            
            % Induction step
            Lred = chol(Mred);
            
            % Get the factor
            O = shiftOperator.zeros(1,m-1);
            N11Sqrt = sqrt(N11);
            L = [N11Sqrt O; N21*inv(N11Sqrt) Lred];
            if (class(L) == "sym")
                L = simplify(L);
            end
        end
        
        function indicies = find(M)
            %FIND Finds first non zero element in a matrix.
            [n,m] = size(M);
            
            indicies = [];
            for i = 1:n
                for j = 1:m
                    op = M(i,j);
                    if ~isZeroOperator(op)
                        indicies(end+1) = i + (j-1)*n;
                        
                    end
                end
            end
        end
        
        function bool = isZeroOperator(m)
            %ISZEROOPERATOR Checks if an operator is the zero operator.
            shiftMatrix = m.shiftCoefficients;
            [n,m] = size(shiftMatrix);
            
            bool = true;
            for i = 1:n
                for j= 1:m
                    if ( shiftMatrix(i,j) ~= 0 )
                        bool = false;
                        return
                    end
                end
            end
            
        end
        
        function bool = onlyZeros(M)
            %ONLYZEROS Checks if a matrix only contains zeros operators.
            %Returns false if any of the elements in the matrix is not a
            %zero operator.
            [n,m] = size(M);
            
            bool = true;
            for i = 1:n
                for j = 1:m
                    op = M(i,j);
                    if ~isZeroOperator(op)
                        bool = false;
                        return
                    end
                end
            end
        end
        
        function bool = eq(A,B)
            %EQ Checks if two matrices of shift operators are equal
            try
                diff = A - B;
            catch
                bool = false;
                return
            end
            bool = onlyZeros(diff);
        end
        
        function bool = ne(A,B)
            %NE Checks if to matrices with shift operators are not equal
            bool = ~(A == B);
        end
        
        function rounded = round(opMatrix,N)
            %ROUND Rounds the shift operator coefficients in a matrix to the nearest decimal or integer.
            %
            %ROUND(A) rounds all elements to the closest integer.
            %
            %ROUND(A,N) rounds all elements in A to N digits right of the
            %decimal point. N has to be positive.
            if nargin == 1
                N = 0;
            end
            [n,m] = size(opMatrix);
            
            rounded = opMatrix;
            for i = 1:n
                for j = 1:m
                    op = opMatrix(i,j);
                    coeff = getCoefficients(op);
                    roundedCoeff = round(coeff, N);
                    roundedOp = shiftOperator(roundedCoeff);
                    rounded(i,j) = roundedOp;
                end
            end
        end
        
        function simplified = simplify(M)
            %SIMPLIFY Simlifies the symbolic expresions in the shift operator coefficients.
            [n,m] = size(M);
            simplified = shiftOperator.zeros(n,m);
            for i = 1:n
                for j = 1:m
                    op = M(i,j);
                    shiftMatrix = getCoefficients(op);
                    if (class(shiftMatrix) == "sym")
                        shiftMatrix = simplify(shiftMatrix);
                    end
                    simplified(i,j) = shiftOperator(shiftMatrix);
                end
            end
        end
        
        function Mpart = monomialCoefficients(M,I,J)
            %MONOMIALCOEFFICIENTS Picks out the coefficients for one type of monomial.
            %MONOMIALCOEFFICIENTS(M,i,j) picks the real valued matrix of
            %coefficients of the  monomial (q')^(j-1)q^(i-1) from M.
            %
            %Example:
            %q = shiftOperator('q');
            %M = [2*q, q'; q+q', 3*q^2];
            %Mq = monomialCoefficients(M,1,2)   %Picks out the matrix [2 0;1 0]
            [n,m] = size(M);
            Mpart = zeros(n,m);
            for i = 1:n
                for j = 1:m
                    % Pick out scalar part of element
                    element = M(i,j);
                    matrixRepresentation = getCoefficients(element);
                    [nShift, mShift] = size(matrixRepresentation);
                    if( I > nShift || J > mShift)
                        Mpart(i,j) = 0;
                    else
                        partOfInterest = matrixRepresentation(I,J);
                        Mpart(i,j) = partOfInterest;
                    end
                end
            end
        end
        
        function d = maxOrder(M)
            %MAXORDER Gives the maximum order of the operators in a matrix.
            %The maximal order of a shift operator is the maximal power of
            %the shift operators making up the operator.
            [n,m] = size(M);
            d = 0;
            for i = 1:n
                for j = 1:m
                    op = M(i,j);
                    opCoeff = getCoefficients(op);
                    opOrder = length(opCoeff) - 1;
                    
                    d = max(d, opOrder);
                end
            end
        end
         
    end
    
    methods (Hidden)
        function operatorSum = plusElement(operator1, operator2)
            %PLUSELEMENT Sums two shiftOperators.
            if( class(operator1) == "double" )
                operator1 = shiftOperator(operator1);
            end
            if (class(operator2) == "double")
                operator2 = shiftOperator(operator2);
            end
            shifts1 = operator1.shiftCoefficients;
            shifts2 = operator2.shiftCoefficients;;
            
            n1 = length(shifts1);
            n2 = length(shifts2);
            if (n2 < n1)
                shifts2(n1,n1) = 0;
            elseif (n1 < n2)
                shifts1(n2,n2) = 0;
            end
            
            operatorSum = shiftOperator(shifts1 + shifts2);
        end
        
        function operatorDifference = minusElement(operator1, operator2)
            %MINUSELEMENT Subtracts one shift operator from another.
            %MINUSELEMENT(op1,op2) performs the calculation op1-op2.
            if( class(operator1) == "double" || class(operator1) == "sym")
                operator1 = shiftOperator(operator1);
            end
            if (class(operator2) == "double"|| class(operator2) == "sym")
                operator2 = shiftOperator(operator2);
            end
            
            shifts1 = operator1.shiftCoefficients;
            shifts2 = operator2.shiftCoefficients;
            
            n1 = length(shifts1);
            n2 = length(shifts2);
            if (n2 < n1)
                shifts2(n1,n1) = 0;
            elseif (n1 < n2)
                shifts1(n2,n2) = 0;
            end
            
            operatorDifference = shiftOperator(shifts1-shifts2);
        end
        
        function operatorProduct = etimes(operator1, operator2)
            %ETIMES Multiplies two shift operators.
            if( class(operator1) == "double" || class(operator1) == "sym")
                operator1 = shiftOperator(operator1);
            end
            if (class(operator2) == "double"|| class(operator2) == "sym")
                operator2 = shiftOperator(operator2);
            end
            shifts1 = operator1.shiftCoefficients;
            shifts2 = operator2.shiftCoefficients;
            
            if (sum(abs(shifts1),'all') == 0 || sum(abs(shifts2),'all') == 0)
                operatorProduct = shiftOperator(0);
                return
            end
            
            [backwardShifts1, forwardShifts1] = size(shifts1);
            [backwardShifts2, forwardShifts2] = size(shifts2);
            
            productShifts = 0;
            for i = 1:backwardShifts1
                for j = 1:forwardShifts1
                    coefficient1 = shifts1(i,j);
                    if (coefficient1 ~= 0)
                        for k = 1:backwardShifts2
                            for l = 1:forwardShifts2
                                coefficient2 = shifts2(k,l);
                                [n,m] = size(productShifts);
                                coeffProduct = coefficient1*coefficient2;
                                if (class(coeffProduct) == "sym")
                                    productShifts = sym(productShifts);
                                end
                                if (j >= k)
                                    if (i > n || j+l-k > m)
                                        productShifts(i,j+l-k) = coeffProduct;
                                    else
                                        productShifts(i,j+l-k) = productShifts(i,j+l-k) + coeffProduct;
                                    end
                                else
                                    if (i+k-j > n || l > m) 
                                        productShifts(i+k-j,l) = coeffProduct;
                                    else
                                        productShifts(i+k-j,l) = productShifts(i+k-j,l) + coeffProduct;
                                    end
                                end
                            end
                        end
                    end
                end
            end
            
            [n,m] = size(productShifts);
            if (n > m)
                productShifts(n,n) = 0;
            elseif (m > n)
                productShifts(m,m) = 0;
            end
            
            operatorProduct = shiftOperator(productShifts);
        end
        
        function product = operatorMatrixProduct(op,M)
            %OPERATORMATRIXPRODUCT Multiplies a matrix with a shift operator from the left.
            %The multiplication is performed element wise
            [n,m] = size(M);
            product = shiftOperator.zeros(n,m);
            for i = 1:n
                for j= 1:m
                    product(i,j) = op*M(i,j);
                end
            end
        end
        
        function product = matrixOperatorProduct(M,op)
            %MATRIXOPERATORPRODUCT Multiplies a matrix with a shift operator from the right
            %The multilpication is performed element wise 
            [n,m] = size(M);
            product = shiftOperator.zeros(n,m);
            for i = 1:n
                for j= 1:m
                    product(i,j) = M(i,j)*op;
                end
            end
        end
        
    end
    
    methods (Static)
        function zOp = zeros(n,m)
            %ZEROS Gives a matrix of zero operators
            %ZERO(n,m) gives a matrix of size n by m
            %ZEROS(n) gives a matrix of size n by n
            
            if nargin == 1
                m = n;
            end
            
            z = shiftOperator(0);
            
            for i = 1:n
                for j = 1:m
                    zOp(i,j) = z;
                end
            end
        end
        
        function I = eye(n)
            %EYE Gives an identity shift operator matrix
            %EYE(n) gives a identity matrix of size n by n
            
            one = shiftOperator(1);
            z = shiftOperator(0);
            I = one;
            for i = 1:n
                for j = i:n
                    if (i == j)
                        I(i,j) = one;
                    else
                        I(i,j) = z;
                    end
                end
            end
        end
        
    end
    
    
end