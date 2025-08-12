using Test
include("ShiftOperator.jl")

@testset "Shift operator class" begin
    @testset "Operators" begin
        @testset "Constructor" begin
            # Build from matrix
            M = [1 2; 0 4];
            try
                op = ShiftOperator(M)
                @test true
            catch
                @test false
            end

            # Build from string
            q = ShiftOperator([0 1])
            qStringF = ShiftOperator("forward")
            @test q == qStringF

            qStringQ = ShiftOperator("q")
            @test q == qStringQ
        end

        @testset "Operator to String" begin
            # Normal
            op = ShiftOperator([0 -3])
            actual = op2str(op)
            expected = "-3q"
            @test actual == expected

            op = ShiftOperator([0 0.1]')
            actual = op2str(op)
            expected = "0.1q'"
            @test actual == expected

            # Zero
            op = ShiftOperator(0)
            actual = op2str(op)
            expected = "0"
            @test actual == expected

            # No 1 coefficients
            op = ShiftOperator([0 0 1])
            actual = op2str(op)
            expected = "q^2"
            @test actual == expected

        end

        @testset "Negative operator" begin
            # Normal operation
            q = ShiftOperator([0 1])
            actual = -q
            expected = ShiftOperator([0 -1])
            @test actual == expected
        end

        @testset "Addition" begin
            # Normal
            op1 = ShiftOperator([0 1])
            op2 = ShiftOperator([0 5])
            actual = op1 + op2
            expected = ShiftOperator([0 6])
            @test actual == expected

            # Multiple terms
            op1 = ShiftOperator([0 1]');
            op2 = ShiftOperator([0 0; 0 3]);
            op3 = ShiftOperator([0 -2]');
            actual = op1 + op2 + op3;
            expected = ShiftOperator([0 0; -1 3])
            @test actual == expected
            
        end

        @testset "Subtraction" begin
            # Normal
            op1 = ShiftOperator([0 1]')
            op2 = ShiftOperator([0 5]')
            op3 = ShiftOperator([0 3])
            actual = op1 - op2 - op3
            expected = ShiftOperator([0 -3; -4 0])
            @test actual == expected
        end

        @testset "Scalar multiplication" begin
            # Normal
            q = ShiftOperator([0 1])
            actual = 2*q
            expected = ShiftOperator([0 2])
            @test actual == expected
            
            # Multiply with zero
            qStar = ShiftOperator([0 1]')
            actual = 0*qStar
            expected = ShiftOperator(0)
            @test actual == expected
            
            # Negative scalar
            actual = -4*q
            expected = ShiftOperator([0 -4])
            @test actual == expected
        end

        @testset "Conjugate" begin
            # Simple
            q = ShiftOperator([0 1])
            qStar = ShiftOperator([0 1]')
            @test qStar == q'
            @test q == qStar'
        end

        @testset "Multipication" begin
            # Basic
            q = ShiftOperator([0 1])
            actual = q*q
            expected = ShiftOperator([0 0 1])
            @test actual == expected

            actual = q*(q')
            expected = ShiftOperator([1])
            @test actual == expected

            actual = (q')*q
            expected = ShiftOperator([0 0; 0 1])
            @test actual == expected
            
            # Normal
            one = ShiftOperator([1])
            actual = (3*q + 5*(q') - (q')*q)*(2*one - q + 2*(q'))
            M = [6 6 -3; 8 -7 1; 10 0 0]
            expected = ShiftOperator(M)
            @test actual == expected
        end

        @testset "Inverse" begin
            # Normal
            q = ShiftOperator([0 1])
            one = ShiftOperator([1])
            op = one + q'*q
            actual = inv(op)
            expected = ShiftOperator([1 0; 0 -0.5])
            @test actual == expected

            # Not in R-infinity
            correctMsg = false
            try
                inv(q)
            catch e
                if e.msg == "operator must be in R-infinity class"
                    correctMsg = true
                end
            end
            @test correctMsg

            # Not invertable
            op = one - q'*q
            correctMsg = false
            try
                inv(op)
            catch e
                if e.msg == "operator is not invertable"
                    correctMsg = true
                end
            end
            @test correctMsg
        end

        @testset "Power" begin
            # Normal
            q = ShiftOperator([0 1])
            actual = q^3
            expected = ShiftOperator([0 0 0 1])
            @test actual == expected

            # Negative
            correctMsg = false
            try
                q^(-1)
            catch e
                if isa(e, MethodError)
                    correctMsg = true
                end
            end
            @test correctMsg

            # Zero
            actual = q^0
            expected = ShiftOperator(1)
            @test actual == expected
            
            # Shift operator
            correctMsg = false
            try
                q^q
            catch e
                if isa(e, MethodError)
                    correctMsg = true
                end
            end
            @test correctMsg
        end

        @testset "Normal definition" begin
            # Normal
            q = ShiftOperator([0 1])
            actual = 4*q' -5*q'*q + 3*q^3
            M = [0 0 0 3; 4 -5 0 0]
            expected = ShiftOperator(M)
            @test actual == expected
        end

        @testset "Square root" begin
            # Normal
            q = ShiftOperator([0 1])
            one = ShiftOperator([1])
            op = 4*one + 5*(q')*q
            actual = sqrt(op)
            expected = 2*one + q'*q
            @test actual == expected

            # Negative
            op = 2*one - 5*(q')*q
            correctMsg = false
            try
                sqrt(op)
            catch e
                if e.msg == "operator must be positive semi-definite"
                    correctMsg = true
                end
            end
            @test correctMsg

            # Not R-infinity
            correctMsg = false
            try
                sqrt(q)
            catch e
                if e.msg == "operator must be in R-infinity class"
                    correctMsg = true
                end
            end
            @test correctMsg
        end
    end

    @testset "Shift operator matrices" begin
        @testset "Zero matrix" begin
            # Two input
            actual = shiftopZeros(ShiftOperator,3,4)
            z = ShiftOperator([0])
            expected = [z z z z; z z z z; z z z z]
            @test actual == expected

            # Two inputs
            actual = shiftopZeros(ShiftOperator, 3)
            expected = [z z z; z z z; z z z]
            @test actual == expected
        end

        @testset "Identity" begin
            # Noraml
            actual = shiftopI(ShiftOperator, 3)
            z = ShiftOperator(0)
            one = ShiftOperator(1)
            expected = [one z z; z one z; z z one]
            @test actual == expected
        end

        @testset "Addition" begin
            # Normal
            q = ShiftOperator([0 1])
            M1 = [4*q' 2*q^2; 2*q' q^2]
            M2 = [2*q q^2; -3*q' q'*q]
            actual = M1 + M2
            expected = [2*q+4*q' 3*q^2; -q' q'*q+q^2]
            @test actual == expected

            # Wrong dimensions
            M3 = [q ; q']
            correctMsg = false
            try
                M1 + M3
            catch e
                if isa(e, DimensionMismatch)
                    correctMsg = true
                end
            end
            @test correctMsg
        end

        @testset "Subtraction" begin
            # Normal
            q = ShiftOperator([0 1])
            M1 = [4*q' 2*q^2; 2*q' q^2]
            M2 = [2*q q^2; -3*q' q'*q]
            actual = M1 - M2
            expected = [-2*q+4*q' q^2; 5*q' -q'*q+q^2]
            @test actual == expected

            # Wrong dimensions
            M3 = [q ; q']
            correctMsg = false
            try
                M1 - M3
            catch e
                if isa(e, DimensionMismatch)
                    correctMsg = true
                end
            end
            @test correctMsg
        end

        @testset "Negative" begin
            # Normal
            I = shiftopI(ShiftOperator, 3)
            actual = -I
            one = ShiftOperator(1)
            z = ShiftOperator(0)
            expected = [-one z z;
                        z -one z;
                        z z -one];
            @test actual == expected
        end

        @testset "Multiplication" begin
            # Normal 2x2 2x2
            q = ShiftOperator([0 1])
            one = ShiftOperator(1);
            M1 = [q' 2*q; one 3*one]
            M2 = [3*q 4*q'; q 2*one]
            actual = M1*M2
            expected = [3*q'*q+2*q^2 4*q+4*(q')^2; 6*q 6*one+4*q']
            @test actual == expected

            # Normal 2x2 2x1
            M3 = [q; q]
            actual = M1*M3
            expected = reshape([q'*q+2*q^2; 4*q], 2, 1)
            @test actual == expected

            # Identity
            I = shiftopI(ShiftOperator,2)
            actual = M1*I
            expected = M1
            @test actual == expected

            # Zero
            O = shiftopZeros(ShiftOperator, 2)
            actual = M1*O
            expected = O
            @test actual == expected
            
            # Wrong dimensions
            correctMsg = false
            try
                M3*M1
            catch e
                if isa(e, DimensionMismatch)
                    correctMsg = true
                end
            end
            @test correctMsg
        end

        @testset "Element wise multiplication" begin
            # Normal
            q = ShiftOperator([0 1])
            z = ShiftOperator(0)
            M = [q 2*q; q^2 q']
            I = shiftopI(ShiftOperator, 2)
            actual = M.*I
            expected = [q z; z q']
            @test actual == expected

            one = ShiftOperator(1)
            N1 = [one; 2*one]
            actual = M.*N1
            expected = [q 2*q; 2*q^2 2*q']
            @test actual == expected

            # Wrong dimensions
            N2 = [q; q; q]
            correctMsg = false
            try
                M.*N2
            catch e
                if isa(e, DimensionMismatch)
                    correctMsg = true
                end
            end
            @test correctMsg
        end

        @testset "Multiplication with scaler" begin
            # Normal
            q = ShiftOperator([0 1])
            M = [q+q' q'; 7*q (q')^2]
            actual = 2*M
            expected = [2*q+2*q' 2*q'; 14*q 2*(q')^2]
            @test actual == expected
        end

        @testset "Multiplication with operator" begin
            # From left
            q = ShiftOperator([0 1])
            one = ShiftOperator(1)
            M = [q+q' q'; 7*q (q')^2]
            actual = q*M
            expected = [one+q^2 one; 7*q^2 q']
            @test actual == expected

            # From right
            actual = M*q
            expected = [q'*q+q^2 q'*q; 7*q^2 (q')^2*q]
            @test actual == expected
        end

        @testset "Adjoint" begin
            # Normal
            q = ShiftOperator([0 1])
            M = [q+q' q'; 7*q (q')^2]
            actual = M'
            expected = [q'+q 7*q'; q q^2]
            @test actual == expected
        end

        @testset "Element wise power" begin
            # Normal
            q = ShiftOperator([0 1])
            M = [q 2*q; q^2 q']
            actual = M.^2
            expected = [q^2 4*q^2; q^4 (q')^2]
            @test actual == expected

            # Operator
            correctMsg = false
            try
                M.^q
            catch e
                if isa(e, MethodError)
                    correctMsg = true
                end
            end
            @test correctMsg

            # Matrix
            correctMsg = false
            A = [1 2; 1 2]
            actual = M.^A
            expected = [q 4*q^2; q^2 (q')^2] 
            @test actual == expected

            #Matrix dimension mismatch
            B = [1 2; 1 2; 1 2]
            correctMsg = false
            try
                M.^B
            catch e
                if isa(e, DimensionMismatch)
                    correctMsg = true
                end
            end
            @test correctMsg
            
            # Negative
            correctMsg = false
            try
                M.^(-1)
            catch e
                if isa(e, MethodError)
                    correctMsg = true
                end
            end
            @test correctMsg
        end

        @testset "Power" begin
            # Normal
            q = ShiftOperator([0 1])
            one = ShiftOperator(1)
            M = [one q; 2*one q']
            actual = M^2;
            expected = [one+2*q q+one; 2*one+2*q' 2*q+(q')^2]
            @test actual == expected

            # Operator
            correctMsg = false
            try
                M^q
            catch e
                if isa(e, MethodError)
                    correctMsg = true
                end
            end
            @test correctMsg

            # Matrix
            correctMsg = false
            A = [1 2; 1 2]
            try
                M^A
            catch e
                if isa(e, MethodError)
                    correctMsg = true
                end
            end
            @test correctMsg
            
            # Negative
            correctMsg = false
            try
                M^(-1)
            catch e
                if isa(e, MethodError)
                    correctMsg = true
                end
            end
            @test correctMsg
        end

        @testset "Max order" begin
            # Normal
            q = ShiftOperator([0 1])
            M = [q^4+q' (q')^5; q q']
            actual = maxOrder(M)
            expected = 5
            @test actual == expected
        end

        @testset "Equalty" begin
            # Normal
            q = ShiftOperator([0 1])
            A = [q' q; 2*q q]
            B = [q' q; q+q q]
            @test A == B

            C = [q' q; q+q 2*q]
            @test A != C

            # Wrong dimnsion
            z = ShiftOperator(0)
            D = [2*q' q; q+q q; z z]
            @test A != D
        end

        @testset "Cholesky" begin
            # String
            q = ShiftOperator([0 1])
            id = ShiftOperator(1)
            z = ShiftOperator(0)
            M = [q' z z; -id q' z; z -id q'; z z -id]
            actual = chol(M)
            expected = [sqrt(2)*id z z;
                        -(1/sqrt(2))*q sqrt(3/2)*id z;
                        z -sqrt(2/3)*q sqrt(4/3)*id]
            diff = actual - expected
            expectedDiff = shiftopZeros(ShiftOperator, 3)
            @test round.(diff) == expectedDiff

            # Not tree
            M = [-id z q'; q' -id z; z q' -id]
            correctMsg = false
            try
                chol(M)
            catch e
                if isa(e, DomainError)
                    correctMsg = true
                end
            end
            @test correctMsg

            # Not elimination ordering
            M = [z z z -id;
                z -id z q';
                z q' -id z;
                -id z q' z;
                q' z z z]
            correctMsg = false
            try
                chol(M)
            catch e
                if isa(e, DomainError)
                    correctMsg = true
                end
            end
            @test correctMsg
        end
    end
end
