classdef signalL2
    %SIGNALL2 
    % This class represents a set of signals in l2+ that can be represented
    % by a tripple such that the signal
    %   w = (C x_0, C A x_0, C A^2 x_0, ...)
    % is represented by the tripple (C, A, x0). The appropriate sizes are
    %   C  - 1 x n (row vector)
    %   A  - n x n (square matrix)
    %   x0 - n x 1 (column vector)
    %
    % The class can be used together with the shiftOperator class. This
    % class can be used when applying such an operator to a signal.
    %
    % Examples:
    % q = shiftOperator('q');
    % op = q^2 + q';
    % syms C A x0;
    % y = signalL2(C, A, x0);
    % yOut = applyOperator(op, y);  % Alterative syntax: op*y
    
    properties
        C       % Output vector
        A       % State matrix
        x0      % Initial condition
    end
    
    methods
        function obj = signalL2(C,A, x0)
            %SIGNALL2 Construct an instance of this class
            %   Detailed explanation goes here
            
            % Check size C
            [nC, mC] = size(C);
            if nC ~= 1
                error("C must be a row vector")
            end
            
            % Check size A
            [nA, mA] = size(A);
            if nA ~= mA
                error("A must be a square matrix")
            elseif mA ~= mC
                error("A and C must have the same number of columns")
            end
            
            % Check size x0
            [nx, mx] = size(x0);
            if mx ~= 1
                error("x0 must be a column vector")
            elseif nx ~= nA
                error("x0 must have the same number of rows as A")
            end
            if class(A) ~= "sym" && class(C) ~= "sym" && class(x0) ~= "sym"
                [C, A, x0] = signalL2.dimensionReduction(C, A, x0);
            end
            obj.C = C;
            obj.A = A;
            obj.x0 = x0;
        end
        
        function C = getC(y)
            %GETC returns the output matrix
            %It gives the vector C of the signal
            %   w = (C x0, C A x0, C A^2 x0, ...)
            C= y.C;
        end
        
        function A = getA(y)
            %GETA returns the state matrix
            %It gives the matrix A of the signal
            %   w = (C x0, C A x0, C A^2 x0, ...)
            A = y.A;
        end
        
        function x0 = getX0(y)
            %GETX0 returns the initial condition
            %It gives the vector x0 of the signal
            %   w = (C x0, C A x0, C A^2 x0, ...)
            x0 = y.x0;
        end
        
        function disp(y)
            %DISP Displays the signals first 3 samples
            [n,m] = size(y);
            if m ~= 1
                error("Signal vectors must be column vectors")
            end
            k = 4;      % Number of samples to show
            if n == 1
                str = firstSamples(y, k);
                disp(str)
            else
                for i = 1:n
                    str = firstSamples(y(i),k);
                    if i == 1
                        str = "[" + str;
                    elseif i == n
                        str = str + "]";
                    end
                    disp(str);
                end
            end
        end
        
        function yPart = firstSamples(y, k)
            %FIRSTSAMPLES Gives the first k samples of signal
            yPart = "(";
            for i = 0:k-1
                sample = y.C*y.A^i*y.x0;
                yPart = yPart + string(sample) + ", ";
            end
            yPart = yPart + "...)";
        end
        
        function signalSum = plus(y1,y2)
            %PLUS Adds two signals
            
            CSum = [y1.C y2.C];
            
            n1 = length(y1.A);
            n2 = length(y2.A);
            ASum = [y1.A zeros(n1,n2); zeros(n2, n1) y2.A];
            
            x0Sum = [y1.x0; y2.x0];
            
            signalSum = signalL2(CSum, ASum, x0Sum);
        end
        
        function signalDiff = minus(y1, y2)
            %MINUS Subtracts two signals
            one = shiftOperator(1);
            signalDiff = y1 + (-one)*y2;
        end
        
        function newY = applyMonomialOperator(cij, i, j, y)
            %APPLYMONOMIALOPERATOR Applies an operator c*(q')^(i-1)*q^(j-1)
            if i < j
                newC = cij*y.C*y.A^(j-i);
                newA = y.A;
                newX0 = y.x0;
            else
                n = length(y.A);
                newC = [zeros(1,(i-1)*n) cij*y.C*y.A^(j-1)];
                newA = zeros(n*i);
                if (class(y.A) == "sym")
                    newA = sym(newA);
                end
                newA(1:n,1:n) = y.A;
                if i > 1
                    I = eye(n*(i-1));
                    newA(n+1:end,1:end-n) = I;
                end
                newX0 = zeros(n*i,1);
                if (class(y.x0) == "sym")
                    newX0 = sym(newX0);
                end
                newX0(1:n,1) = y.x0;
            end
            newY = signalL2(newC, newA, newX0);
        end
        
    end
    
    methods (Static)
        function [Cred, Ared, x0red] = dimensionReduction(C, A, x0)
            % Reduces the dimension to the minimum necessary dimension
            % for observability
            Wo = obsv(A, C);
            n = rank(Wo);
            if n == length(A)
                Cred = C;
                Ared = A;
                x0red = x0;
                return
            end
            if n == 0
                Cred = 0;
                Ared = 0;
                x0red = 0;
                return
            end
            
            [~,~,v] = svd(Wo);
            Ared = v'*A*v;
            Ared = Ared(1:n, 1:n);
            Cred = C*v;
            Cred = Cred(1, 1:n);
            x0red = v'*x0;
            x0red = x0red(1:n, 1);
        end
    end
end

