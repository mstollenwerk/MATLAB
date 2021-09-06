function [ eparam, tstats, logL, optimoutput ] = matvFRiesz2est( data_mat, x0, varargin )
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
[p,~,N] = size(data_mat);
p_ = p*(p+1)/2;
%%
if ~isempty(varargin) && strcmp(varargin{1},'EstimationStrategy:MethodOfMoments')
    %% Optimization
    
    mean_data_mat = mean(data_mat,3);
    U = cholU(mean_data_mat);
    
    obj_fun = @(df) matvFRiesz2like( U/matvFRiesz2expmat(df(1:p),df(p+1:p+p))*U', df(1:p), df(p+1:p+p), data_mat );

    if isempty(x0)  
        x0 = [ 2*p.*ones(p,1); 2*p.*ones(p,1) ];
    end

    lb = [(p-(1:p)+2), (0:p-1)]'; %conjectured

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
        'df_1', tstats(1 : p), ...               
        'df_2', tstats(p + 1 : p + p), ...               
        'all', [ NaN(p_,1); tstats ] ...
    );
    %% nLogL, logLcontr and eparam
    [nLogL, logLcontr, ~, ~, eparam ] = obj_fun( eparam );

    aic = 2*nLogL + 2*(p_ + p + p);
    bic = 2*nLogL + log(N)*(p_ + p + p);
    logL = struct(...
        'logL', -nLogL, ...
        'logLcontr', logLcontr, ...
        'bic', bic, ...
        'aic', aic ...
    );
else
    %% Optimization
    obj_fun = @(param) matvFRiesz2like( [], [], [], data_mat, param );

    if isempty(x0)  
        x0 = [ vech(chol(mean(data_mat,3)/2/p, 'lower')); 2*p.*ones(p,1); 2*p.*ones(p,1) ];
    end

    lb = [-inf(p_,1); (p-(1:p)+2)'; (0:p-1)'];

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
        'df_1', tstats(p_ + 1 : p_ + p), ...               
        'df_2', tstats(p_ + p + 1 : p_ + p + p), ...               
        'all', tstats ...
    );
    %% nLogL, logLcontr and eparam
    [nLogL, logLcontr, ~, ~, eparam ] = obj_fun( eparam );

    aic = 2*nLogL + 2*(p_ + p + p);
    bic = 2*nLogL + log(N)*(p_ + p + p);
    logL = struct(...
        'logL', -nLogL, ...
        'logLcontr', logLcontr, ...
        'bic', bic, ...
        'aic', aic ...
    );
end

end
