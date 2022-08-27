function [ eparam, nLogL, h, tau, V_m ] = mf2_garch_estim(dta)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
m=63; % 21, 63, 126, 252

obj_fun = @(param) mf2_garch_likeRec(param, dta);

% x0 = [0.007; 0.85; 0.01; 0.07; 2; 10];
x0 = [0.007; 0.85; 0.01; 0.07; 1; 10];

% lb = [0, 0, 0, 0, 1, 1];
lb = [0, 0, 0, 0, 0, 1];
% ub = [1, 1, inf, 1, inf, inf];
ub = [1, 1, inf, 1, 1, inf];

options = optimoptions(@fmincon, "Display", "iter-detailed", 'MaxFunEval', 1e3);
eparam = fmincon( ...
    obj_fun, x0, ...
    [], [], [], [], lb, ub, [], ...
    options ...
);

[ nLogL, h, tau, V_m ] = mf2_garch_likeRec(eparam, dta);

end

