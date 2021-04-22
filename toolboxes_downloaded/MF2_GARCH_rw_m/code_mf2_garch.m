
% Code for estimating an (MF)^2 GARCH-rw-m Model

% open the return data
load('sp500_long.mat'); 
y = RET_SPX_LONG;

% parameters: alpha, gamma, beta, lambda_0, lambda_1, lambda_2
% start values
%param_init = [0.02; 0.007; 0.14 ; 0.85; 0.01; 0.05; 0.93 ]; % m = 63 
param_init = [0.02; 0.007; 0.14 ; 0.85; 0.01; 0.07; 0.91 ]; % m = 63


% constraints:
% alpha >=0
% alpha + gamma/2 + beta <=1
% lambda_1 >= 0
% lambda_1 + lambda_2 <=1
A = [  ...
    0.0 -1.0  0.0  0.0  0.0  0.0  0.0; ...
    0.0 1.0  0.5  1.0  0.0  0.0  0.0; ...
    0.0 0.0  0.0  0.0  0.0  -1.0  0.0; ...
    0.0 0.0  0.0  0.0  0.0  1.0  1.0...
];
b = [  0.0; 1.0; 0.0; 1.0 ];

Aeq = [];
beq = [];
nlc = [];
% bounds:
LB = [-1; 0.0 ; -0.5; 0.0; 0.000001; 0.0; 0.0 ];
UB = [ 1; 1.0 ; 0.5; 1.0; 0.2; 1.0; 1.0 ];
% options = optimset('fmincon');
% options = optimoptions(@fmincon,'Algorithm','sqp') ;
options = optimoptions(@fmincon, "Algorithm", "active-set", "Display", "iter-detailed");

% choice for m
m=63; % 21, 63, 126, 252

% estimation
[ coeff, ll, exitFlag, ~, ~, ~, hessian ] = fmincon('likelihood_mf2_garch', param_init, A, b, Aeq, beq, LB, UB, nlc, options, y, m)


[ e, h, tau, V_m ] = mf2_garch_core(coeff, y, m);

