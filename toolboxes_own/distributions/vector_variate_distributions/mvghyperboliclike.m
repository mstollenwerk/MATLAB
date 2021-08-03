function [ nLogL, logLcontr, varargout ] = ...
    mvghyperboliclike( ...
        mu_, Sigma_, beta_, alpha_, lambda_, delta_, data_mat, varargin ...
    )
%
%
% USAGE:
%   
%
% INPUTS:
%       mu_    - \xi in [1]
%       Sigma_ - V in [1]
%
% OUTPUTS:
%   
%  See also 
%
% COMMENTS:
%   
% REFERENCES:
%      [1] Blaesild and Jensen (1981) - Multivariate Distributions
%          of Hyperbolic Type
%      [2] Kotz, S., Balakrishnan, N., & Johnson, N. L. (2004). 
%          Continuous multivariate distributions, Volume 1: 
%          Models and applications. John Wiley & Sons.
%          p. 212.
%      [3] Vanduffel and Yao (2017) - A stein type lemma for the 
%          multivariate generalized hyperbolic distribution.
%      [4] McNeil, A. J. , Frey, R. , & Embrechts, P. (2005). 
%          Quantitative risk management . Princeton University Press,
%          Princeton, NJ.
%
% Michael Stollenwerk
% michael.stollenwerk@live.com
% 20.09.2020
%
% DEPENDENCIES:
%

warning('This file is not finished. Paramter checks and boundary cases need to be implemented. See ref [1].')

narginchk(7,8);
nargoutchk(0,5);

[N,p] = size(data_mat);
p_ = p*(p+1)/2;
%% Parameters
if nargin == 8
    all_param = varargin{1};
    mu_ = all_param(1:p);
    Sigma_ = ivechchol(all_param( p+1 : p+p_ ));
    beta_ = all_param(p + p_ + 1 : p + p_ + p);
    alpha_ = all_param(p + p_ + p + 1);
    lambda_ = all_param(p + p_ + p + 2);
    delta_ = all_param(p + p_ + p + 3);
end

if nargout == 5
    param.mu_ = mu_;
    param.Sigma_ = Sigma_;
    param.beta_ = beta_;    
    param.alpha_ = alpha_;
    param.lambda_ = lambda_;
    param.delta_ = delta_;
    param.all = [mu_
                 vech(chol(Sigma_,'lower'))
                 beta_ 
                 alpha_ 
                 lambda_ 
                 delta_]';
    
    varargout{3} = param;
end
%% Log-Likelihood
logLcontr = NaN(N,1);

eta_ = sqrt(alpha_^2 - beta_*Sigma_*beta_);

log_norm_const = ...
    lambda_ * log(eta_/delta_) ...
    + (p/2 - lambda_)*log(alpha_) ...
    - p/2*log(2*pi) ...
    - log(besselk(lambda_,delta_*eta_));

invSig = inv(Sigma_);
parfor ii = 1:N

    x = data_mat(ii,:) - mu_;
    q = delta_^2 + x*invSig*x';
    
    log_kernel = ...
        log(besselk( lambda_ - p/2, alpha_*sqrt(q) )) ...
        - (p - 2*lambda_)/4 * log(q) ...
        + beta_*x';
    
    logLcontr(ii) = log_norm_const + log_kernel;
    
%     if isinf(logLcontr(ii))
%         nLogL = logLcontr(ii);
%         return
%     end
%     
%     if isnan(logLcontr(ii))
%         nLogL = logLcontr(ii);
%         return
%     end
    
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
