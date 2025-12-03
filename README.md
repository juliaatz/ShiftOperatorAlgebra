# Shift operator algebra
This is the code associated with the PAPER.
It implements an algebra for shift operators that is used to derive
factorised version of the Linear Quadratic Regulator (LQR).
The code is available in two languages; MATLAB and Julia.

Each subdirectory contains the following (language specific) files
- **Examples** containes the code for the examples in the paper.
- **Factor_lqr** is a function that derives a factorised LQR algorithimically.
- **Factor_lqr_formula** is a function that computes a factorised LQR
    by the formula in Lemma X from the paper.
- **M2ss** is a function that derives a state space formulation from a
    matrix of shift operators (when certain conditions are fullfilled).
- **PlotGraph** is a function that plots the graph associated with a matrix of shift operators
- **RandomTree** is a function that creates a matrix of shift operators
    representing a random tree shaped graph.
- **ShiftOperator** is the class that defines the shift operator object
- **ShiftOperatorTutorial** is a texted base tutorial that explaines how
    the shift operator class can be used.
- **TimeVaryingNecessary** shows that a time varying factorisation is needed
    in accordance with lemma X.

For installation and more language specific instructions, check the READMEs in the
corresponding subdirectories.
