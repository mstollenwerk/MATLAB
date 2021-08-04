function [ eparam, tstats, logL, optimoutput ] = wishest( data_mat, varargin )
%WISHEST Estimates the parameters of the Wishart distribution.
%
% USAGE:
%   [ eparam, tstats, logL, optimoutput ] = wishest( data_mat, varargin )
%
% INPUTS:
%   DATA_MAT     - p by p by N array of covariance matrix data.
%   VARARGIN     - 1. x0 - Optimization Starting Value.
%                - 2.    - Option-value pairs for optimizer.
%
% OUTPUTS:
%   EPARAM       - Struct of estimated parameters with following fields:
%       SIGMA_    - p by p by N array of scale matrices.
%       DF        - Degree of freedom paramter.
%   TSTATS       - Struct of t-stats with the fields as eparam.
%   LOGL         - Struct of estimated parameters with following fields:
%       LOGL      - Log-likelihood value.
%       LOGLCONTR - Log-likelihood contributions. 
%   OPTIMOUTPUT  - Optimization Output from fmincon with fields:
%                 - estimation_time 
%
% COMMENTS:
%   This is not W(df, Sigma/df)!!!
%
%  See also WISHLIKE WISHQLIKE WISHPDF WISHVARMAT WISHVARMAT_OLD
%
% REFERENCES:
%      [1] Gupta, Nagar (1999) - Matrix Variate Distributions, p. 87.

% Michael Stollenwerk
% michael.stollenwerk@live.com
% 30.01.2020
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
    if isfield(x0,'Sigma_')
        Sigma_ = x0.Sigma_;
    else
        Sigma_ = eye(p);
    end    
    if isfield(x0,'df')
        df = x0.df;
    else
        df = 2*p;
    end 
    x0 = [ vech(chol(Sigma_, 'lower')); df ]'; % k*(k+1)/2 + 1
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
    'MaxFunctionEvaluations', 1e5, ...
    'UseParallel', 'always' ...
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
[nLogL, logLcontr, Sigma_, df ] = like_wrapper(eparam, data_mat, p, N);

eparam = struct( ...
    'Sigma_', Sigma_, ...
    'df', df, ...
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
    'df', tstats(p*(p+1)/2 + 1) ...
);
%% Likelihood wrapper function
function [nLogL, logLcontr, Sigma_, df ] = like_wrapper(param, data_mat, p, N) 

p_ = p*(p+1)/2;

Sigma_ = ivech(param(1:p_), 'lower');
Sigma_ = Sigma_ * Sigma_';
df = param(p_+1);

[ nLogL, logLcontr ] = wishlike( repmat(Sigma_,1,1,N), df, data_mat );

end
end

