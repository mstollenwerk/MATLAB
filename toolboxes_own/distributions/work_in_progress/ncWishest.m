function [ eparam, tstats, logL, optimoutput ] = ncWishest( data_mat, varargin )
%WISHEST Estimates the parameters of the Wishart distribution.
%
% USAGE:
%   [ eparam, tstats, logL, optimoutput ] = ncWishest( data_mat, varargin )
%
% INPUTS:
%   DATA_MAT     - p by p by N array of covariance matrix data.
%   VARARGIN     - 1. x0 - Optimization Starting Value.
%                - 2. MHG_PRECISION - Number of Jack Functions used to 
%                                     approximate the matrix valued 
%                                     hypergeometric function. This input 
%                                     is inversly related to computing time.
%                - 3.    - Option-value pairs for optimizer.
%
% OUTPUTS:
%   EPARAM       - Struct of estimated parameters with following fields:
%       DF            - Vector of length N. Degrees of freedom parameters 
%                       of non-central Wishart distribution.
%       SIGMA_        - p by p by N array of parameter matrices, 
%                       regulates the covariance. 
%       OMEGA_        - p by p by N array of non-centrality matrices.
%   TSTATS       - Struct of t-stats with the fields as eparam.
%   LOGL         - Struct of estimated parameters with following fields:
%       LOGL          - Log-likelihood value.
%       LOGLCONTR     - Log-likelihood contributions. 
%       BIC
%       AIC
%   OPTIMOUTPUT  - Optimization Output from fmincon with fields:
%                 - estimation_time 
%
% COMMENTS:
%   OMEGA_ is NOT equal to the THETA matrix in Gupta, Nagar (1999), but 
%   THETA_ = inv(SIGMA_)*OMEGA_.
%
% REFERENCES:
%   [1] Gupta, Nagar (1999) - Matrix Variate Distributions, p. 114.
%
% DEPENDENCIES:
%
%  See also NCWISHLIKE
%
% Michael Stollenwerk
% michael.stollenwerk@live.com
% 30.01.2020
%% Error checking
narginchk(1,inf);
if nargin==1
    x0 = [];
    mhg_precision = [];
elseif nargin == 2
    x0 = varargin{1};
    mhg_precision = [];
else
    x0 = varargin{1};
    mhg_precision = varargin{2};
end
[p,p2,N] = size(data_mat);
if p~=p2
    error('Covariance matrix data must be square.')
end
%% Optimization
% Objective Function ------------------------------------------------------
fun = @( param ) like_wrapper(param, data_mat, p, N, mhg_precision);
% x0 ----------------------------------------------------------------------
if isempty(x0) || isstruct(x0)
    if isfield(x0,'Sigma_')
        Sigma_ = x0.Sigma_;
    else
        Sigma_ = eye(p);
    end    
    if isfield(x0,'Omega_')
        Omega_ = x0.Omega_;
    else
        Omega_ = eye(p);
    end        
    if isfield(x0,'df')
        df = x0.df;
    else
        df = 2*p;
    end 
    x0 = [ vech(chol(Sigma_, 'lower')); vech(chol(Omega_, 'lower')); df ]'; % k*(k+1) + 1
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
options = optimoptions(options, varargin{3:end});
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
[nLogL, logLcontr, Sigma_, Omega_, df ] = like_wrapper(eparam, data_mat, p, N, mhg_precision);

eparam = struct( ...
    'Sigma_', Sigma_, ...
    'Omega_', Omega_, ...
    'df', df ...
);

aic = 2*nLogL + 2*(numel(x0));
bic = 2*nLogL + log(p*(p+1)/2*N)*(numel(x0));
logL = struct(...
    'logL', -nLogL, ...
    'logLcontr', logLcontr, ...
    'bic', bic, ...
    'aic', aic ...
);
tstats = struct(...
    'Sigma_', tstats(1:p*(p+1)/2),...
    'Omega_', tstats(p*(p+1)/2+1:p*(p+1)),...    
    'df', tstats(p*(p+1) + 1) ...
);
%% Likelihood wrapper function
function [nLogL, logLcontr, Sigma_, Omega_, df ] = like_wrapper(param, data_mat, p, N, mhg_precision) 

p_ = p*(p+1)/2;

Sigma_ = ivech(param(1:p_), 'lower');
Sigma_ = Sigma_ * Sigma_';
Omega_ = ivech(param(p_+1:p_+p_), 'lower');
Omega_ = Omega_ * Omega_';
df = param(p_+p_+1);

[ nLogL, logLcontr ] = ncWishlike( data_mat, ones(N,1)*df, repmat(Sigma_,1,1,N), repmat(Omega_,1,1,N), mhg_precision);

end
end

