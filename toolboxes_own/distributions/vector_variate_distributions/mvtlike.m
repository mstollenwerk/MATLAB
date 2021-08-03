function [ nLogL, logLcontr, varargout ] = ...
    mvtlike( mu_, Sigma_, nu_, data_mat, varargin )
%MVTLIKE Negative log-likelihood and score of the multivariate t-distribution.
%
% USAGE:
%   [ nLogL, logLcontr, varargout ] = mvtlike( mu_, Sigma_, nu_, data_mat, varargin )
%
% INPUTS:
%   DATA_MAT- Array (p by 1). Data.
%   MU_     - Array (p by 1). Mean parameters.
%   SIGMA_  - Array (p by p). Symmetric p.d. parameter matrix, 
%             regulates covariance. 
%   NU_     - Double. Degrees of freedom parameter.
%
% OUTPUTS:
%   NLOGL   - Double. Negative log-likelihood value.
%   SCORE   - Struct. Fields as parameter names. Contains derivatives
%             of log-likelihood (matvFlike) w.r.t. parameters.
%
% See also MVTSCORE MVTRND
%
% COMMENTS:
%
%   Derivatives were derived using matrixcalculus.org
%
% REFERENCES:
%   [1] Azzalini, Adelchi. The skew-normal and related families. 
%           Vol. 3. Cambridge University Press, 2013.
%   [2] My distributions tex. This is the distribuion NOT dividing by df-2!
%
% DEPENDENCIES:
%  gammaln
%
% Michael Stollenwerk
% michael.stollenwerk@live.com
% 01.06.2020

narginchk(4,5);
nargoutchk(0,5);

if size(data_mat,2) ==1
    data_mat = data_mat';
end

[N,p] = size(data_mat);
p_ = p*(p+1)/2;
%% Parameters
if nargin == 5
    all_param = varargin{1};
    if size(all_param,1) < size(all_param,2)
        all_param = all_param';
    end
    mu_ = all_param(1:p);
    Sigma_ = ivechchol(all_param( p+1 : p+p_ ));
    nu_ = all_param(p + p_ + 1);
end

if isempty(mu_)
    mu_ = zeros(p,1);
end

if nargout == 5
    param.mu_ = mu_;
    param.Sigma_ = Sigma_;
    param.lambda_ = nu_;
    param.all = [mu_
                 vech(chol(Sigma_,'lower'))
                 nu_]';
    
    varargout{3} = param;
end
%%
if nu_ <= 1
    warning('Degree of freedom parameter is smaller than 1. For these, Mean and Covariance Matrix do not exist')
elseif nu_ <= 2
    warning('Degree of freedom parameter is smaller than 2. For these, Covariance Matrix does not exist')
end
%% Log-Likelihood
logLcontr = NaN(N,1);

log_norm_const = ...
    gammaln( (nu_ + p)/2 ) ...
    - gammaln(nu_/2) ...
    - p/2*log( nu_*pi ) ...
    - log(det(Sigma_))/2;

invSigma = inv(Sigma_);
for ii = 1:N
    
    x = data_mat(ii,:) - mu_';
    
    log_kernel = ...
        -(nu_+p)/2*log( 1 + x*invSigma*x'/nu_ );
    
    logLcontr(ii) = log_norm_const + log_kernel;
    
end
nLogL = -sum(logLcontr);
%% Scores
if nargout >= 3
    
    scoreSig = NaN(N,p_);
    score_nu = NaN(N,1);
    
    for ii = 1:N
        
        x = data_mat(ii,:) - mu_';
       
        S = 1/2*(nu_ + p)/(nu_ + x*invSigma*x')*invSigma*x'*x*invSigma ...
            - 1/2*invSigma ;
        
        % Apply vec operator.
        S = S(:);
        
        % Accounting for symmetry of Sigma_:
        D = Dmatrix(p);
        scoreSig(ii,:) = D'*S;

        score_nu(ii) = NaN;
    
    end
    
    score.Sigma_ = scoreSig;
    score.nu_ = score_nu;
    
    varargout{1} = score;
    
end
%% Fisher Info
if nargout == 4
    
    varargout{2} = NaN;
    
end
end
