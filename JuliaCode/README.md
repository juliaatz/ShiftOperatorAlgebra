# Julia README

## Installation
Installed the code by clonning the git repository by running
    gitclone https://github.com/juliaatz/ShiftOperatorAlgebra.git
in the directory you want the code to be stored.

To Download Julia check the instructions [here](https://julialang.org/downloads/).

## Using the code
To run the code, open a terminal and navigate to the cloned repository and then the JuliaCode
subdirectory. After that start the Julia REPL by typing `julia`.

In order to have the same package versions as was originally intended a PKg environoment is avaiable.
To activate it enter the Pkg mode in the Julia REPL by pressing `]`. There type the following commands
    activate .
    instantiate
To exit the Pkg mode press `Ctrl + c`.

In order to be able to use a code file write `include("fileName.jl")` to be able to acess the code from the REPL.
As an example
  include("ShiftOperators.jl")
allows you to use the ShiftOperator class. To get started with it it is recommended to check out the tutorial `ShiftOperatorTutorial.md`.
To run the code from the examples in the paper, write
  include("examples.jl")
