function [eparam, tstats, logL, optimoutput] = matvgWisht3est(X, x0, df_n, varargin)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

%% Input Checking
% Will be added later
[p,~,N] = size(X);
p_ = p*(p+1)/2;
df_n_ = df_n*(df_n+1)/2;
%% Optimization
obj_fun = @(param) matvgWisht3like([],[],[],X,param);

if isempty(x0)  
    x0 = [ vechchol(mean(X,3)*(5-2)/5); vech(eye(df_n)); 5 ]';
end

% Restrictions-------------------------------------------------------------
lb = [ -inf(p_ + df_n_,1); 1];
% Optimization-------------------------------------------------------------
[eparam,optimoutput] = ...
    my_fmincon( ...
        obj_fun, ...
        x0,[],[],[],[],lb,[],[],varargin{:} ...
    );

%% tstats
%[VCV,A,B,scores,hess,gross_scores] = robustvcv(fun, eparam, 3);
[VCV,scores,gross_scores] = vcv(obj_fun, eparam);

tstats = eparam./sqrt(diag(VCV));

tstats = struct(...
    'vechcholSigma_', tstats(1:p_), ...
    'vechcholTheta_', tstats(p_ + 1 : p_ + df_n_), ...
    'df_t', tstats(p_ + df_n_ + 1), ...           
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

