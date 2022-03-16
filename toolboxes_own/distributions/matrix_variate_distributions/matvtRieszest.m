function [ eparam, tstats, logL, optimoutput ] = matvtRieszest( R, x0, varargin )
%BESSELWISHEST
%
% USAGE:
%   
%
% INPUTS:
%   
%
% OUTPUTS:
%   
%  See also 
%
% COMMENTS:
% 
% REFERENCES:
%      [1] Stollenwerk (2020)
%
% Michael Stollenwerk
% michael.stollenwerk@live.com
% 11.02.2021
%
% DEPENDENCIES:
%
narginchk(2,inf);
[p,~,N] = size(R);
p_ = p*(p+1)/2;
if ~isempty(varargin) && strcmp(varargin{1},'EstimationStrategy:MethodOfMoments')
    %% Optimization
    
    if isempty(x0)  
        x0 = [ 2*p.*ones(p,1); 2*p ];
    end

    lb = [ (0:p-1)'; 2 ];

    avgR = mean(R,3);
    obj_fun = @(df,perm_) matvstRieszlike( avgR(perm_,perm_,:), ... 
                                           df(1:p), df(p+1), ...
                                           R(perm_,perm_,:) );
    
    [eparam,optimoutput] = ...
        fmincon_Rieszperm(...
            p, ...
            obj_fun, ...
            x0, ...
            [],[],[],[],lb,[],[],varargin{2:end} ...
        );
    
    obj_fun_opt_perm = @(df) obj_fun(df,optimoutput.perm_);
    %% tstats
    %[VCV,A,B,scores,hess,gross_scores] = robustvcv(fun, eparam, 3);
    [VCV,scores,gross_scores] = ...
        vcv(obj_fun_opt_perm, eparam);

    tstats = eparam./sqrt(diag(VCV));

    tstats = struct(...
        'Sigma_', NaN(p_,1), ...
        'n', tstats(1:p), ...
        'nu', tstats(p+1), ...               
        'all', [NaN(p_,1); tstats] ...
    );
    %% nLogL, logLcontr and eparam
    [nLogL, logLcontr, ~, ~, eparam ] = obj_fun_opt_perm( eparam );
    eparam.perm_ = optimoutput.perm_;

    aic = 2*nLogL + 2*(p_ + p + 1);
    bic = 2*nLogL + log(N)*(p_ + p + 1);
    logL = struct(...
        'logL', -nLogL, ...
        'logLcontr', logLcontr, ...
        'bic', bic, ...
        'aic', aic ...
    );
else
    %% Optimization
    warning("Optimization over all parameters is done without optimization over asset ordering!")    
    obj_fun = @(param) matvtRieszlike( [], [], [], R, param );

    if isempty(x0)  
        x0 = [ vech(chol(mean(R,3)/2/p, 'lower')); 2*p.*ones(p,1); 2*p ];
    end

    lb = [-inf(p_,1); (0:p-1)'; 2];

    [eparam,optimoutput] = ...
        my_fmincon(...
            obj_fun, ...
            x0, ...
            [],[],[],[],lb,[],[],varargin{:} ...
        );
    %% tstats
    %[VCV,A,B,scores,hess,gross_scores] = robustvcv(fun, eparam, 3);
    [VCV,scores,gross_scores] = ...
        vcv(obj_fun, eparam);

    tstats = eparam./sqrt(diag(VCV));

    tstats = struct(...
        'Sigma_', tstats(1:p_), ...
        'n', tstats(p_ + 1 : p_ + p), ...               
        'nu', tstats(p_ + p + 1), ... 
        'all', tstats ...
    );
    %% nLogL, logLcontr and eparam
    [nLogL, logLcontr, ~, ~, eparam ] = obj_fun( eparam );

    aic = 2*nLogL + 2*(p_ + p + 1);
    bic = 2*nLogL + log(N)*(p_ + p + 1);
    logL = struct(...
        'logL', -nLogL, ...
        'logLcontr', logLcontr, ...
        'bic', bic, ...
        'aic', aic ...
    );
end

end