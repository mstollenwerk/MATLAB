function [eparam, tstats, logL, optimoutput] = matvncsFBOPSest(X, mhg_precision, rank, x0, varargin)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

%% Input Checking
% Will be added later
[p,~,N] = size(X);
p_ = p*(p+1)/2;
%% Optimization
% Objective function
obj_fun = ...
    @(param) matvncFBOPSlike( ...
        ivechchol(param(1:p_)), ...
        reshape(param(p_+1:p_+p*rank),p,rank) ...
            *reshape(param(p_+1:p_+p*rank),p,rank)', ...
        param(p_+p*rank + 1), ...
        param(p_+p*rank + 2), ...
        mhg_precision, X ...
    );
% x0
if isempty(x0)
    %"Pre-Optimization" based on central F.
    x0 = matvFest(X, []);
    %"Pre-Optimization" to get good "Sigma_ to Omega_ ratio".
    meanX = x0.Sigma_*x0.df_1/(x0.df_2 - 2);
    weight = .1:.1:1;
    y = NaN(length(weight),1);
    for ii = 1:length(weight) 
        y(ii) = matvncFBOPSlike( ...
            meanX*weight(ii), ...
            meanX*(1 - weight(ii))*x0.df_1, ...
            x0.df_1, x0.df_2 + p - 1, mhg_precision, X ...
        ); 
    end
    
    x0 = [vechchol(meanX*weight(y==min(y)));
          ones(p*rank,1)*sqrt(mean(mean(meanX))*(1-weight(y==min(y)))*x0.df_1);
          x0.df_1;
          x0.df_2 + p - 1];
      
% You should start far away from Omega_ = 0-Matrix!
%     x0 = [vechchol(eye(p).*.01);vechchol(mean(X,3)*2*p);2*p;2*p];
%     x0 = [vechchol(mean(X,3));vechchol(eye(p).*.1);2*p;2*p];
end
% Minimize objective function
[eparam,optimoutput] = my_fminunc( obj_fun, x0, varargin{:} );
%% tstats
%[VCV,A,B,scores,hess,gross_scores] = robustvcv(fun, eparam, 3);
[VCV,scores,gross_scores] = vcv(obj_fun, eparam);

tstats = eparam./sqrt(diag(VCV));

tstats = struct(...
    'Sigma_', tstats(1:p_), ...
    'Omega', tstats(p_+1 : p_+p*rank), ...
    'df_1', tstats(p_+p*rank + 1), ...           
    'df_2', tstats(p_+p*rank + 2), ...     
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
