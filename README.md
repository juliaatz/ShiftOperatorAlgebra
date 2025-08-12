# shift-operator-algebra
This is the code associated with the PAPER.
It implements an algebra for shift operators that is used to derive
factorised version of the Linear Quadratic Regulator (LQR).
The code is available in two languages; MATLAB and Julia.

Each subdirectory contains the following (language specific) files
- **Examples** containes the code for the examples in the paper.
- **Factor_lqr** is a function that derives a factorised LQR algorithimically.
- **Factor_lqr_formula** is a function that derives a factorised LQR
    by the formula in Lemma X from the paper.
- **M2ss** is a function that derives a state space formulation from a
    matrix of shift operators (when certain conditions are fullfilled).
- **PlotGraph** is a function that plots the graph associated with a matrix of shift operators
- **RandomTree** is a function that creates a matrix of shift operators
    representing a random tree shaped graph.
- **ShiftOperator** is the class that defines the shift operator object
- **ShiftOperatorTutorial** is a texted base tutorial that explaines how
    the shift operator class can be used.

The MATLAB code has been developed in MATLAB 2020b. For the Julia code
a Pkg envirnoment is available through the `Project.toml` and
`Manifest.toml` files. These can be activated Pkg mode in the Julia REPL
(press "]" after writing julia in the terminal) with the following commands:

        activate JuliaCode
        instantiate
