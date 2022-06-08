function [ eparam, tstats, logL, optimoutput ] = matviWishest( data_mat, x0, varargin )
%MATVIWISHEST
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
%      [1]
%
% DEPENDENCIES:
%
% Michael Stollenwerk
% michael.stollenwerk@live.com

narginchk(2,inf);
[p,~,N] = size(data_mat);
p_ = p*(p+1)/2;
if ~isempty(varargin) && strcmp(varargin{1},'EstimationStrategy:MethodOfMoments')
    %% Optimization
    mean_data_mat = mean(data_mat,3);
    obj_fun = @(df) matviWishlike( mean_data_mat*(df-p-1), df, data_mat );
    
    if isempty(x0)  
        x0 = p+10;
    end

    lb = p+1;

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
        'df', tstats, ...               
        'all', [NaN(p_,1);tstats] ...
    );
    %% nLogL, logLcontr and eparam
    [nLogL, logLcontr, ~, ~, eparam ] = obj_fun( eparam );

    aic = 2*nLogL + 2*(p_+1);
    bic = 2*nLogL + log(N)*(p_+1);
    logL = struct(...
        'logL', -nLogL, ...
        'logLcontr', logLcontr, ...
        'bic', bic, ...
        'aic', aic ...
    );
else
    %% Optimization
    obj_fun = @(param) matviWishlike( [], [], data_mat, param );
    if isempty(x0)  
        x0 = [ vech(chol(mean(data_mat,3)*(2*p - p - 1), 'lower')); 2*p ]';
    end

    lb = [-inf(p_,1);p+1];

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
        'df', tstats(p_ + 1), ...               
        'all', tstats ...
    );
    %% nLogL, logLcontr and eparam
    [nLogL, logLcontr, ~, ~, eparam ] = obj_fun( eparam );

    aic = 2*nLogL + 2*numel(x0);
    bic = 2*nLogL + log(N)*numel(x0);
    logL = struct(...
        'logL', -nLogL, ...
        'logLcontr', logLcontr, ...
        'bic', bic, ...
        'aic', aic ...
    );
end

end
