function [ nLogL, logLcontr, varargout ] = ...
    mvsymhyperboliclike( ...
        mu_, Sigma_, alpha_, delta_, data_mat, varargin ...
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
%      [1] Blaesild and Jensen (1981) - Multivariate Distributions
%          of Hyperbolic Type. Eq. 9.
%      [2] McNeil, A. J. , Frey, R. , & Embrechts, P. (2005). 
%          Quantitative risk management . Princeton University Press,
%          Princeton, NJ.
%
% Michael Stollenwerk
% michael.stollenwerk@live.com
% 22.09.2020
%
% DEPENDENCIES:
%
error('This is still wrong. Use Quantitative Risk Management Book to correct.)
narginchk(5,6);
nargoutchk(0,5);

[N,p] = size(data_mat);
p_ = p*(p+1)/2;
%% Parameters
if nargin == 6
    all_param = varargin{1};
    mu_ = all_param(1:p);
    Sigma_ = ivechchol(all_param( p+1 : p+p_ ));
    alpha_ = all_param(p + p_ + 1);
    delta_ = all_param(p + p_ + 2);
end

if nargout == 5
    param.mu_ = mu_;
    param.Sigma_ = Sigma_;
    param.alpha_ = alpha_;
    param.delta_ = delta_;
    param.all = [mu_
                 vech(chol(Sigma_,'lower'))
                 alpha_ 
                 delta_]';
    
    varargout{3} = param;
end
%% Log-Likelihood
logLcontr = NaN(N,1);

log_norm_const = ...
    (p+1)/2 * log(alpha_/delta_) ...
    - (p-1)/2*log(2*pi) ...    
    - log(2*alpha_) ...
    - log(besselk( (p+1)/2, delta_*alpha_ ));

invSig = inv(Sigma_);
for ii = 1:N

    x = data_mat(ii,:) - mu_;
    
    log_kernel = -alpha_*sqrt(delta_^2 + x*invSig*x');
    
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
