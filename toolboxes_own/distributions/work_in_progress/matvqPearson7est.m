function [eparam, tstats, logL, optimoutput] = matvqPearson7est(X, x0, varargin)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

%% Input Checking
% Will be added later
[p,~,N] = size(X);
p_ = p*(p+1)/2;
%% Optimization
obj_fun = @(param) matvqPearson7like([],[],[],[],X,param);

if isempty(x0)
    x0 = [vech(chol(eye(p),'lower'));4*p;p^2+p;2*p];
end
% Restrictions-------------------------------------------------------------
A = [zeros(1,p_) p/2 -1 0];
b = 0;
lb = [ -inf(p_,1); 2*p; -inf; 0];
% Optimization-------------------------------------------------------------
    [eparam,optimoutput] = ...
    my_fmincon( ...
        obj_fun, ...
        x0,A,b,[],[],lb,[],[],varargin{:} ...
    );

%% tstats
%[VCV,A,B,scores,hess,gross_scores] = robustvcv(fun, eparam, 3);
[VCV,scores,gross_scores] = vcv(obj_fun, eparam);

tstats = eparam./sqrt(diag(VCV));

tstats = struct(...
    'Sigma_', tstats(1:p_), ...
    'n', tstats(p_ + 1), ...           
    'q', tstats(p_ + 2), ...     
    'r', tstats(p_ + 3), ...  
    'all', tstats ...
);
%% nLogL, logLcontr and eparam  
[nLogL, logLcontr, ~, ~, eparam] = obj_fun(eparam);

aic = 2*nLogL + 2*numel(x0);
bic = 2*nLogL + log(N)*numel(x0);
logL = struct(...
    'logL', -nLogL, ...
    'logLcontr', logLcontr, ...
    'bic', bic, ...
    'aic', aic ...
);
end

