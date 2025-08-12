%% Set up
clear; close all
q = shiftOperator('q');
one = shiftOperator(1);
z = shiftOperator(0);
r = 0.5;                % Discount factor

%% Introduction example
% Matrix representing the transportation routs and delays
m = 4;
p = m+1;
M = shiftOperator.zeros(p,m);
for i = 1:m
    M(i,i) = q'*r;
    M(i+1,i) = -one;
end

% plotGraph(M);     % Uncomment if you want to plot the graph


% Get the state space representation
sys = M2ss(M);
A = sys.A;
C = sys.C;
B = sys.B;

% Find control law
[K1, K2] = factor_lqr_formula(A, C, M, r)

% Check validity
K = K1\K2;
Q = C'*C;
R = zeros(m);
[~, K_idare, ~] = idare(A*r,B,Q,R);
Kerror = K - K_idare

%% Simple example
% Matrix representing the transportation routs and delays
M = [q'*r z z;
    -one q'*r z;
    z -one q'*r;
    z z -one];

% plotGraph(M);     % Uncomment if you want to plot the graph

[p,m] = size(M);

% Get the state space representation
sys = M2ss(M);
A = sys.A;
C = sys.C;
B = sys.B;

% Find control law
[K1, K2] = factor_lqr_formula(A, C, M, r)

% Check validity
K = K1\K2;
Q = C'*C;
R = zeros(m);
[~, K_idare, ~] = idare(A*r,B,Q,R);
Kerror = K - K_idare

%% Tree 1
% Matrix representing the transportation routs and delays
M = [q'*r z z z;
    z q'*r z z;
    -one z q'*r z;
    z -one -one q'*r;
    z z z -one];

% plotGraph(M);     % Uncomment if you want to plot the graph

[p,m] = size(M);

% Get the state space representation
sys = M2ss(M);
A = sys.A;
C = sys.C;
B = sys.B;

% Find control law
[K1, K2] = factor_lqr_formula(A, C, M, r)

% Check validity
K = K1\K2;
Q = C'*C;
R = zeros(m);
[~, K_idare, ~] = idare(A*r,B,Q,R);
Kerror = K - K_idare

%% Tree 2
% Matrix representing the transportation routs and delays
M = [q'*r z z z;
    z -one z z;
    -one z -one z;
    z q'*r q'*r q'*r;
    z z z -one];

% plotGraph(M);     % Uncomment if you want to plot the graph

[p,m] = size(M);

% Get the state space representation
sys = M2ss(M);
A = sys.A;
C = sys.C;
B = sys.B;

% Find control law
[K1, K2] = factor_lqr_formula(A, C, M, r)

% Check validity
K = K1\K2;
Q = C'*C;
R = zeros(m);
[~, K_idare, ~] = idare(A*r,B,Q,R);
Kerror = K - K_idare

%% String
% Matrix representing the transportation routs and delays
M = [-one z z z;
    q'*r -one z z;
    z q'*r q'*r z;
    z z -one q'*r;
    z z z -one];

% plotGraph(M);     % Uncomment if you want to plot the graph

[p,m] = size(M);

% Get the state space representation
sys = M2ss(M);
A = sys.A;
C = sys.C;
B = sys.B;

% Find control law
[K1, K2] = factor_lqr(A, C, M, r)

% Check validity
K = K1\K2;
Q = C'*C;
R = zeros(m);
[~, K_idare, ~] = idare(A*r,B,Q,R);
Kerror = K - K_idare

%% Tree 3
rng(2)

% Creat random matrix representing the transportation routs and delays
m = 20;
mixedDirection = false;
outflow = -one;
inflow = q'*r;
M = randomTree(m, mixedDirection, outflow, inflow);

% plotGraph(M);     % Uncomment if you want to plot the graph

[p,m] = size(M);


% Get the state space representation
sys = M2ss(M);
A = sys.A;
C = sys.C;
B = sys.B;

% Find control law
[K1, K2] = factor_lqr(A, C, M, r)

% Check validity
K = K1\K2;
Q = C'*C;
R = zeros(m);
[~, K_idare, ~] = idare(A*r,B,Q,R);
Kerror = K - K_idare

