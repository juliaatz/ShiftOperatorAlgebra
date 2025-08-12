# About this tutorial

This is a tutorial for how to use the class shiftOperator in Julia.
This is not meant as an exhaustive explanation of the class but just to get
you started. For more details on the shiftOperator class, see the documentation.
To see how the class can be used to derive a factorised linear quadratic
regulator, see examples.jl.

# Basic shift operators

All shift operators consist of two building blocks: the forward and the
backward shift. They are denoted `q` and `q'`.

## Construction 

Before starting, make sure that the ShiftOperator class is loaded with
`include("ShiftOperator.jl")`.

The simplest way to define the forward shift operator is

    q = ShiftOperator("q")

All shift operators are represented by a matrix containing the multiplication
table of its coefficients. You can also define a shift operator through
this matrix, e.g

    q = ShiftOperator([0 1])        # Gives the forward shift
    qAdj = ShiftOperator([0 1]')    # Gives the backward shift
    op = ShiftOperator([1 2; 3 4])  # Gives the operator 1 + 2q + 3q' + 4q'q

## Arithmetic

Ones you have created a shift operator object, you can quickly make more by
combining them with arithmetic operations.

### Basic arithmetic

All the basic arithmetic operations are available with operator overloading:
addition, subtraction, multiplication with scalar, multiplication with
other shift operators, power. You can also compute the adjoint with the ' operator.

    q = ShiftOperator("q")
    op = 2*q + q'*q + 4*q^2

### Extra arithmetic for R-infinity

The class R-infinity is defined and explained properly in the PAPER.
It is the set of shift operators that can be written as 
$y = \sum_{i = 0}^N a_i (q')^i q^i$.

For operators in R-infinity we also provide the square rote operation and
the inverse (if they exists).

    q = ShiftOperator("q")
    id = ShiftOperator(1)
    op = id + 3*q'*q
    sqrtOp = sqrt(op)
    invOp = inv(op)

# Matrices of shift operators

Shift operators can be put into matrices, e.g.

    q = ShiftOperator("q")
    M = [2*q q'; 5*q 3*q'*q]

You can also create some useful matrices with the `shiftopZeros` and
`shiftopI` functions,
e.g.

    O2 = shiftopZeros(ShiftOperator,2)
    O23 = shiftopZeros(ShiftOperator,2,3)
    I3 = shiftopI(ShiftOperator,3)

## Arithmetic

For matrices of shift operators the basic arithmetic operations are 
available with operator overloading; addition, subtraction, multiplication,
element wise multiplication, multiplication with scalar, power,
element wise power, adjoint.

    q = ShiftOperator("q")
    id = ShiftOperator(1)
    A = [id q; 2*q 3*q']
    B = [4*id q; 3*q q]
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

    q = ShiftOperator('q')
    id = ShiftOperator(1)
    z = ShiftOperator(0)
    M = [q'*r z z;
        -id q'*r z;
        z -id q'*r;
        z z -id]
    L = chol(M)         # L*L' = M'*M
