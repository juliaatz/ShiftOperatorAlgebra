include("ShiftOperator.jl")
include("factor_lqr_formula.jl")
include("factor_lqr.jl")
include("M2ss.jl")
include("randomTree.jl")
using DelimitedFiles                # Required to load the large tree example
using Graphs, Plots, GraphRecipes   # Required if you wan to plot the graphs
include("plotGraph.jl")             # Required if you wan to plot the graphs

q = ShiftOperator("q")
id = ShiftOperator(1)
z = ShiftOperator(0)
r = 0.5
outflow_pattern = -id
inflow_pattern = q'*r

#------------------- Simple example -------------------------
# Matrix representing the transportation routs and delays
M = [q'*r z z;
     -id q'*r z;
     z -id q'*r;
     z z -id]
#plotGraph(M, outflow_pattern, inflow_pattern)   # Uncomment to plot the graph

p,m = size(M)
A,B,C,D = M2ss(M)

K1, K2 = factor_lqr_formula(A,C,M,r)
display(K1)
print("\n")
display(K2)
print("\n")
print("\n")

#------------------- Tree 1 -------------------------
# Matrix representing the transportation routs and delays
M = [q'*r z z z;
    z q'*r z z;
    -id z q'*r z;
    z -id -id q'*r;
    z z z -id]
#plotGraph(M, outflow_pattern, inflow_pattern)   # Uncomment to plot the graph

p,m = size(M)
A,B,C,D = M2ss(M)

K1, K2 = factor_lqr_formula(A,C,M,r)
display(K1)
print("\n")
display(K2)
print("\n")
print("\n")

#------------------- Tree 2 -------------------------
# Matrix representing the transportation routs and delays
M = [q'*r z z z;
    z -id z z;
    -id z -id z;
    z q'*r q'*r q'*r;
    z z z -id]
#plotGraph(M, outflow_pattern, inflow_pattern)   # Uncomment to plot the graph

p,m = size(M)
A,B,C,D = M2ss(M)

K1, K2 = factor_lqr_formula(A,C,M,r)
display(K1)
print("\n")
display(K2)
print("\n")
print("\n")

#------------------- String -------------------------
# Matrix representing the transportation routs and delays
M = [-id z z z;
    q'*r -id z z;
    z q'*r q'*r z;
    z z -id q'*r;
    z z z -id]
#plotGraph(M, outflow_pattern, inflow_pattern)   # Uncomment to plot the graph

p,m = size(M)
A,B,C,D = M2ss(M)

K1, K2 = factor_lqr(A,C,M,r)
display(K1)
print("\n")
display(K2)
print("\n")
print("\n")

#------------------- Big  tree -------------------------
# The big tree in the artikel was generated with the randomTree function
# in Matlab. We provide an extra example bellow that uses Julias randomTree
# function.
M0 = readdlm("M0-big-tree.txt",',','\n')
M1 = readdlm("M1-big-tree.txt",',','\n')
M = M0*id + M1*q'
#plotGraph(M, outflow_pattern, inflow_pattern)   # Uncomment to plot the graph

p,m = size(M)
A,B,C,D = M2ss(M)

K1, K2 = factor_lqr(A,C,M,r)
display(K1)
print("\n")
display(K2)
print("\n")
print("\n")

#------------------- Random  tree -------------------------
m = 10
mixedDirections = false

M = randomTree(m, mixedDirections, outflow_pattern, inflow_pattern)
#plotGraph(M, outflow_pattern, inflow_pattern)   # Uncomment to plot the graph

p,m = size(M)
A,B,C,D = M2ss(M)

K1, K2 = factor_lqr(A,C,M,r)
display(K1)
print("\n")
display(K2)
print("\n")
print("\n")
