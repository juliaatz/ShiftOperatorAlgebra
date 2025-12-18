# About this tutorial

This is a tutorial for how to use the class shiftOperator as well as the
class signalL2 in Matlab.
This is not meant as an exhaustive explanation of the class but just to get
you started. For more details on the shiftOperator class, see the documentation.
To see how the class can be used to derive a factorised linear quadratic
regulator, see examples.m.

# Basic shift operators

All shift operators consist of two building blocks: the forward and the
backward shift. They are denoted `q` and `q'`.

## Construction 

Before starting, make sure that the shiftOperator class is in the current
Matlab path.

The simplest way to define the forward shift operator is

    q = shiftOperator("q")

All shift operators are represented by a matrix containing the multiplication
table of its coefficients. You can also define a shift operator through
this matrix, e.g

    q = shiftOperator([0 1])        % Gives the forward shift
    qAdj = shiftOperator([0; 1])    % Gives the backward shift
    op = shiftOperator([1 2; 3 4])  % Gives the operator 1 + 2q + 3q' + 4q'q

## Arithmetic

Ones you have created a shift operator object, you can quickly make more by
combining them with arithmetic operations.

### Basic arithmetic

All the basic arithmetic operations are available with operator overloading:
addition, subtraction, multiplication with scalar, multiplication with
other shift operators, power. You can also compute the adjoint with the ' operator.

    q = shiftOperator("q")
    op = 2*q + q'*q + 4*q^2

### Extra arithmetic for R-infinity

The class R-infinity is defined and explained properly in the PAPER.
It is the set of shift operators that can be written as 
$y = \sum_{i = 0}^N a_i (q')^i q^i$.

You can check if an operator is in the set with the function `isRinfOperator`.

For operators in R-infinity we also provide the square rote operation and
the inverse (if they exists).

    q = shiftOperator("q");
    op = 1 + 3*q'*q;
    sqrtOp = sqrt(op)
    invOp = inv(op)

# Matrices of shift operators

Shift operators can be put into matrices, e.g.

    q = shiftOperator("q");
    M = [2*q q'; 5*q 3*q'*q]

You can also create some useful matrices with the `zeros` and `eye` functions,
e.g.

    O2 = shiftOperator.zeros(2)
    O23 = shiftOperator.zeros(2,3)
    I3 = shiftOperator.eye(3)

## Arithmetic

For matrices of shift operators the basic arithmetic operations are 
available with operator overloading; addition, subtraction, multiplication,
element wise multiplication, multiplication with scalar, power,
element wise power, adjoint.

    q = shiftOperator("q");
    one = shiftOperator(1);
    A = [one q; 2*q 3*q'];
    B = [4*one q; 3*q q];
    sumAB = A + B
    diff = A - B
    mult = A*B
    scalMult = 2*A
    elMult = A.*B
    pow = A^2
    elPow = A.^2
    adjA = A'

Besides basic arithmetic operations, the Cholesky algorithm from the PAPER
is implemented in the function `chol`. Note that the matrix needs to have
the specific structure required from the PAPER.

    q = shiftOperator('q')
    one = shiftOperator(1)
    z = shiftOperator(0)
    M = [q'*r z z;
        -one q'*r z;
        z -one q'*r;
        z z -one]
    L = chol(M)         % L*L' = M'*M

## Applying operators

To apply the operator you need a signal of the class signalL2, also provided
in this repository. The signal class represents signals that can be represented
by a tripple $(C,A,x_0)$ where the signal is
        $y = (C x_0, C A x_0, C A^2 x_0, C A^3 x_0, ...)$.

The following example shows how an operator can be applied to signal;

    q = shiftOperator('q');
    op = q^2 + q';
    y = signalL2([1 0], [2 1; 0 1], [0;1]);
    yOut = applyOperator(op, y);  % Alterative syntax: op*y
