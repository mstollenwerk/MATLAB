function [ eparam, tstats, logL, optimoutput ] = mFest( data_mat, varargin )
%MFEST Estimates the parameters of the matrix variate F distribution.
%
% USAGE:
%   [ eparam, tstats, logL, optimoutput ] = mFest( data_mat, varargin )
%
% INPUTS:
%   DATA_MAT     - p by p by N array of covariance matrix data.
%   VARARGIN     - 1. x0 - Optimization Starting Value.
%                - 2.    - Option-value pairs for optimizer.
%
% OUTPUTS:
%   EPARAM       - Struct of estimated parameters with following fields:
%      DF_1          - First degrees of freedom parameter of rescaled matrix variate F distribution.
%      DF_2          - Second degrees of freedom parameter of rescaled matrix variate F distribution.
%      SIGMA_        - p by p parameter matrix, regulates the covariance. 
%   TSTATS       - Struct of t-stats with the fields as eparam.
%   LOGL         - Struct of estimated parameters with following fields:
%       LOGL      - Log-likelihood value.
%       LOGLCONTR - Log-likelihood contributions. 
%   OPTIMOUTPUT  - Optimization Output from fmincon with fields:
%                 - estimation_time 
%
%  See also MFLIKE MFRND
%
% COMMENTS:
%   This is the scaled matrix variate F distribution. Scaled s.th.
%   expecatation equals Sigma_.
%
% REFERENCES:
%      [1] 
%
% Michael Stollenwerk
% michael.stollenwerk@live.com
% 04.02.2020
%
% DEPENDENCIES:
%   MVGAMMALN MFLIKE
%% Error checking
narginchk(1,inf);
if nargin==1
    x0 = [];
else
    x0 = varargin{1};
end
[p,p2,N] = size(data_mat);
if p~=p2
    error('Covariance matrix data must be square.')
end
%% Optimization
% Objective Function ------------------------------------------------------
fun = @( param ) like_wrapper(param, data_mat, p, N);
% x0 ----------------------------------------------------------------------
if isempty(x0) || isstruct(x0) 
    if isfield(x0,'df_1')
        df_1 = x0.df_1;
    else
        df_1 = 2*p;
    end
    if isfield(x0,'df_2')
        df_2 = x0.df_2;
    else
        df_2 = 2*p;
    end  
    if isfield(x0,'Sigma_')
        Sigma_ = x0.Sigma_;
    else
        Sigma_ = mean(data_mat,3)*(df_1 - p - 1)/df_2;
    end      
    x0 = [ vech(chol(Sigma_, 'lower')); df_1; df_2 ]'; % k*(k+1)/2 + 2
    try
        fun(x0);
    catch
        error('x0 is not valid.')
    end
else
    error('x0 must be empty or struct.')
end
% Restrictions-------------------------------------------------------------
% Optimoptions-------------------------------------------------------------
options = optimoptions('fminunc',...
    'MaxFunctionEvaluations', 1e5 ...
);
options = optimoptions(options, varargin{2:end});
% Optimization-------------------------------------------------------------
warning('off','MATLAB:nearlySingularMatrix') %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic
[eparam,~,exitflag,optimoutput,grad,hessian] = fminunc(...
    fun,x0,options);
optimoutput.estimation_time = toc;
optimoutput.exitflag = exitflag;
optimoutput.grad = grad;
optimoutput.hessian = hessian;
warning('on','MATLAB:nearlySingularMatrix') %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Tstats-------------------------------------------------------------------
%[VCV,A,B,scores,hess,gross_scores] = robustvcv(fun, eparam, 3);
[VCV,scores,gross_scores] = vcv(fun, eparam);
tstats = eparam./sqrt(diag(VCV)');
% Output Creation----------------------------------------------------------
[nLogL, logLcontr, Sigma_, df_1, df_2 ] = like_wrapper(eparam, data_mat, p, N);

eparam = struct( ...
    'Sigma_', Sigma_, ...
    'df_1', df_1, ...
    'df_2', df_2, ...    
    'all', eparam ...
);

aic = 2*nLogL + 2*(numel(x0));
bic = 2*nLogL + log(N)*(numel(x0));
logL = struct(...
    'logL', -nLogL, ...
    'logLcontr', logLcontr, ...
    'bic', bic, ...
    'aic', aic ...
);
tstats = struct(...
    'Sigma_', tstats(1:p*(p+1)/2),...
    'df_1', tstats(p*(p+1)/2 + 1), ...
    'df_2', tstats(p*(p+1)/2 + 2) ...    
);
%% Likelihood wrapper function
function [nLogL, logLcontr, Sigma_, df_1, df_2 ] = like_wrapper(param, data_mat, p, N) 

p_ = p*(p+1)/2;

Sigma_ = ivech(param(1:p_), 'lower');
Sigma_ = Sigma_ * Sigma_';
df_1 = param(p_+1);
df_2 = param(p_+2);

[ nLogL, logLcontr ] = mFlike( ...
    data_mat, ones(1,N)*df_1, ones(1,N)*df_2, repmat(Sigma_,1,1,N) ...
);

end
end

