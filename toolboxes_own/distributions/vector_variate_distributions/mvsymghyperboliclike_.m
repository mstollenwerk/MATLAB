function [ nLogL, logLcontr, varargout ] = ...
    mvsymghyperboliclike( ...
        mu_, Sigma_, lambda_, chi_, psi_, data_mat, varargin ...
    )
%
%
% USAGE:
%   
%
% INPUTS:
%
% OUTPUTS:
%   
%  See also 
%
% COMMENTS:
%   
% REFERENCES:
%      [1] McNeil, A. J. , Frey, R. , & Embrechts, P. (2005). 
%          Quantitative risk management . Princeton University Press,
%          Princeton, NJ.
%      [2] Blaesild and Jensen (1981) - Multivariate Distributions
%          of Hyperbolic Type
%      [3] Kotz, S., Balakrishnan, N., & Johnson, N. L. (2004). 
%          Continuous multivariate distributions, Volume 1: 
%          Models and applications. John Wiley & Sons.
%          p. 212.
%      [4] Vanduffel and Yao (2017) - A stein type lemma for the 
%          multivariate generalized hyperbolic distribution.
%
% Michael Stollenwerk
% michael.stollenwerk@live.com
% 09.10.2020
%
% DEPENDENCIES:
%

narginchk(6,7);
nargoutchk(0,5);

[N,p] = size(data_mat);
p_ = p*(p+1)/2;
%% Parameters
if nargin == 7
    all_param = varargin{1};
    mu_ = all_param(1:p);
    Sigma_ = ivechchol(all_param( p+1 : p+p_ ));
    lambda_ = all_param(p + p_ + 1);
    chi_ = all_param(p + p_ + 2);
    psi_ = all_param(p + p_ + 3);
end

if nargout == 5
    param.mu_ = mu_;
    param.Sigma_ = Sigma_;
    param.lambda_ = lambda_;
    param.chi_ = chi_;
    param.delta_ = psi_;
    param.all = [mu_
                 vech(chol(Sigma_,'lower'))
                 lambda_ 
                 chi_ 
                 psi_]';
    
    varargout{3} = param;
end
%% Log-Likelihood
logLcontr = NaN(N,1);

eta_ = sqrt(alpha_^2 - beta_*Sigma_*beta_);

log_norm_const = ...
    - lambda_*log(sqrt(chi_*psi_)) ...
    + p/2*log(psi_) ...
    - p/2*log(2*pi) ...
    - log(besselk(lambda_,sqrt(chi_*psi_)))     ...
    - .5*log(det(Sigma_));

invSig = inv(Sigma_);
parfor ii = 1:N

    x = data_mat(ii,:) - mu_;
    sqrt_q = sqrt( (chi_ + x'*invSig*x)*psi_ );
    
    log_kernel = ...
        log(besselk( lambda_ - p/2, sqrt_q )) ...
        - (p/2 - lambda_) * log(sqrt_q);
    
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
