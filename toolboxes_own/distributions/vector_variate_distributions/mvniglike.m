function [ nLogL, logLcontr, varargout ] = ...
    mvniglike( mu_, Sigma_, alpha_, beta_, data_mat, varargin )
%
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
% Michael Stollenwerk
% michael.stollenwerk@live.com
% 19.09.2020
%
% DEPENDENCIES:
%

narginchk(5,6);
nargoutchk(0,5);

[N,p] = size(data_mat);
p_ = p*(p+1)/2;

warning('This distribution seems like just a reparameterized t-distribution.')
%% Parameters
if nargin == 6
    all_param = varargin{1};
    mu_ = all_param(1:p);
    Sigma_ = ivechchol(all_param( p+1 : p+p_ ));
    alpha_ = all_param(p + p_ + 1);
    beta_ = all_param(p + p_ + 2);
end

if nargout == 5
    param.mu_ = mu_;
    param.Sigma_ = Sigma_;
    param.alpha_ = alpha_;
    param.beta_ = beta_;
    param.all = [mu_, vech(chol(Sigma_,'lower'))', alpha_, beta_];
    
    varargout{3} = param;
end
%% Log-Likelihood
logLcontr = NaN(N,1);

log_norm_const = ...
    alpha_ * log(beta_) ...
    + gammaln( alpha_ + p/2 ) ...
    - gammaln( alpha_ ) ...
    - p/2*log(2*pi) ...
    - .5*log(det(Sigma_));

invSig = inv(Sigma_);
for ii = 1:N
    
    x = data_mat(ii,:) - mu_;

    log_kernel = -(alpha_ + p/2)*log(x*invSig*x' / 2 + beta_);
    
    logLcontr(ii) = log_norm_const + log_kernel;
    
end
nLogL = -sum(logLcontr);
%% Scores
if nargout >= 3
    
    varargout{1} = NaN;
    
end
%% Fisher Info
if nargout == 4
    
    varargout{2} = NaN;
    
end
end
