function [ eparam, innov, condVar, nLogL, lls, fcst_condVar ] = garch_x_est( returndata, X, X_g, P, Q, x0 )
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

% Michael Stollenwerk
% michael.stollenwerk@live.com
% 03.05.2017
% 20.06.2017: Changed to GlobalSearch algo to robustly find optimum.

%% Likelihood Optimization
lb = [-inf*ones(size(X,2),1); -ones(1+P+Q+size(X_g,2),1)*.1];
ub = [inf*ones(size(X,2),1); ones(1+P+Q+size(X_g,2),1)];
if isempty(x0)
    x0 = [.1*ones(1,size(X,2)) .01 .1*ones(1,P)/P .85*ones(1,Q)/Q zeros(1,size(X_g,2))];
end

opts = optimoptions('fmincon','Display','iter-detailed','MaxFunctionEvaluations',1e5,'StepTolerance',1e-5,'Display','off','Algorithm','sqp'); % ,'UseParallel','always'
obj_fcn = @(param)garch_x_like(param, returndata, X, X_g, 1, 1);
problem = createOptimProblem('fmincon','x0',x0,...
    'objective',obj_fcn,'lb',lb,'ub',ub,'options',opts);
gs = GlobalSearch;
[eparam,~] = run(gs,problem);
% lb = [-inf*ones(size(X,2),1); -ones(1+P+Q+size(X_g,2),1)];
% ub = [inf*ones(size(X,2),1); ones(1+P+Q+size(X_g,2),1)];
% warning('off', 'optimlib:checkbounds:PadLbWithMinusInf')
% try % If provided x0 is invalid starting point for solver use default x0
%     eparam = fmincon(@(param)garch_x_like(param, returndata, X, X_g, 1, 1), x0,[],[],[],[],lb,[],[],options);
% catch
%     eparam = fmincon(@(param)garch_x_like(param, returndata, X, X_g, 1, 1), x0,[],[],[],[],lb,[],[],options);
% end
[nLogL, lls, innov, condVar, fcst_condVar] = garch_x_like(eparam, returndata, X, X_g, 1, 1);
end