function [eparam, tstats, logL, optimoutput] = ...
    mvsymhyperbolicest(data_mat, x0, varargin)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

error('This is still wrong. Use Quantitative Risk Management Book to correct.)
narginchk(2,inf);
[N,p] = size(data_mat);
p_ = p*(p+1)/2;
%% Optimization
if isempty(x0)
    x0 = [zeros(p,1);
          vech(chol(cov(data_mat)/2,'lower'));
          1;
          1]';
end

lb = [ -inf(p+p_,1); 0; 0 ]';

[eparam,optimoutput] = ...
    my_fmincon( ...
        @(param) mvsymhyperboliclike( [], [], [], [], data_mat, param ), ...
        x0, ...
        [],[],[],[],lb,[],[],varargin{:} ...
    );
%% tstats
%[VCV,A,B,scores,hess,gross_scores] = robustvcv(fun, eparam, 3);
[VCV,scores,gross_scores] = ...
    vcv(@(param) mvsymhyperboliclike( [], [], [], [], data_mat, param ), eparam);

tstats = eparam./sqrt(diag(VCV));

tstats = struct(...
    'mu_', tstats(1:p), ...
    'Sigma_', tstats( p+1 : p+p_ ), ...
    'alpha_', tstats(p + p_ + 1), ...
    'delta_', tstats(p + p_ + 3), ...    
    'all', tstats ...
);
%% nLogL, logLcontr and eparam
[nLogL, logLcontr, ~, ~, eparam ] = ...
    mvsymhyperboliclike( [], [], [], [], data_mat, eparam );

aic = 2*nLogL + 2*numel(x0);
bic = 2*nLogL + log(N)*numel(x0);
logL = struct(...
    'logL', -nLogL, ...
    'logLcontr', logLcontr, ...
    'bic', bic, ...
    'aic', aic ...
);
end
