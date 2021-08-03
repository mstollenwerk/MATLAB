function [ eparam, tstats, logL, optimoutput ] = matvtWishest( X, x0, varargin )
%MATVTWISHEST
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
%   This is the distribution of the quadratic form XX', where the columns
%   of the p by df_n matrix X follow multivariate central t distributions
%   with the same degree of freedom nu and the same dispersion matrix 
%   Sigma_.
%   
% REFERENCES:
%      [1]                    
%
% Michael Stollenwerk
% michael.stollenwerk@live.com
% 11.08.2020
%
% DEPENDENCIES:

narginchk(2,inf);
[p,~,N] = size(X);
p_ = p*(p+1)/2;
%%
if ~isempty(varargin) && strcmp(varargin{1},'EstimationStrategy:MethodOfMoments')
    %% Optimization
    
    mean_X = mean(X,3);
    obj_fun = @(df) matvtWishlike( mean_X*(df(2)-2)/df(2)/df(1), df(1), df(2), X );
    
    if isempty(x0)  
        x0 = [ 2*p; 5 ];
    end

    lb = [p-1;2];

    [eparam,optimoutput] = ...
        my_fmincon(...
            obj_fun, ...
            x0, ...
            [],[],[],[],lb,[],[],varargin{2:end} ...
        );
    %% tstats
    %[VCV,A,B,scores,hess,gross_scores] = robustvcv(fun, eparam, 3);
    [VCV,scores,gross_scores] = vcv(obj_fun, eparam);

    tstats = eparam./sqrt(diag(VCV));

    tstats = struct(...
        'Sigma_', NaN(p_,1), ...
        'df_n', tstats(1), ...           
        'df_t', tstats(2), ...     
        'all', [NaN(p_,1); tstats] ...
    );
    %% nLogL, logLcontr and eparam
    [nLogL, logLcontr, ~, ~, eparam ] = obj_fun( eparam );

    aic = 2*nLogL + 2*(p_+2);
    bic = 2*nLogL + log(N)*(p_+2);
    logL = struct(...
        'logL', -nLogL, ...
        'logLcontr', logLcontr, ...
        'bic', bic, ...
        'aic', aic ...
    );    
else
    %% Optimization
    obj_fun = @(param) matvtWishlike( [], [], [], X, param );

    if isempty(x0)  
        x0 = [ vech(chol(mean(X,3)*(5-2)/5/2/p, 'lower')); 2*p; 5 ];
    end

    lb = [-inf(p_,1);p-1;2];

    [eparam,optimoutput] = ...
        my_fmincon(...
            obj_fun, ...
            x0, ...
            [],[],[],[],lb,[],[],varargin{:} ...
        );
    %% tstats
    %[VCV,A,B,scores,hess,gross_scores] = robustvcv(fun, eparam, 3);
    [VCV,scores,gross_scores] = vcv(obj_fun, eparam);

    tstats = eparam./sqrt(diag(VCV));

    tstats = struct(...
        'Sigma_', tstats(1:p_), ...
        'df_n', tstats(p_ + 1), ...           
        'df_t', tstats(p_ + 2), ...     
        'all', tstats ...
    );
    %% nLogL, logLcontr and eparam
    [nLogL, logLcontr, ~, ~, eparam ] = obj_fun( eparam );

    aic = 2*nLogL + 2*(p_+2);
    bic = 2*nLogL + log(N)*(p_+2);
    logL = struct(...
        'logL', -nLogL, ...
        'logLcontr', logLcontr, ...
        'bic', bic, ...
        'aic', aic ...
    );
end

end
