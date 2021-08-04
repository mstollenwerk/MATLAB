function [ nLogL, logLcontr, varargout ] = ...
    mvsymghyperboliclike( ...
        mu_, Sigma_, lambda_, chi_, psi_, data_mat, varargin ...
    )
%
%
% USAGE:    [ nLogL, logLcontr, varargout ] = ...
%               mvsymghyperboliclike( mu_, Sigma_, lambda_, chi_, psi_, data_mat, varargin )
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
%      [1] McNeil, A. J. , Frey, R. , & Embrechts, P. (2015). 
%          Quantitative risk management . Princeton University Press,
%          Princeton, NJ, p.186.
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
% 20.09.2020
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
    if size(all_param,1) < size(all_param,2)
        all_param = all_param';
    end
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
    param.psi_ = psi_;
    param.all = [mu_
                 vech(chol(Sigma_,'lower'))
                 lambda_ 
                 chi_ 
                 psi_]';
    
    varargout{3} = param;
end
%% Boundary Cases
if psi <= eps && psi >= -eps && ...
   lambda_ <= chi_ + eps && lambda_ >= chi - eps
    [ nLogL, logLcontr, score, fisherInfo ] = ...
        mvtlike(mu_, Sigma_, chi_, data_mat);
    varargout{1} = score;
    varargout{2} = fisherInfo;
    return
end
%More to come. See e.g. McNeil, Frey, Embrechts, p. 186.
%% Log-Likelihood
logLcontr = NaN(N,1);

log_norm_const = ...
    p/2 * log(psi_/2/pi) ...
    - lambda_/2 * log(chi_*psi_) ...
    - log(besselk( lambda_, sqrt(chi_*psi_) )) ...
    - .5*log(det(Sigma_));

invSig = inv(Sigma_);
parfor ii = 1:N

    x = data_mat(ii,:) - mu_';
    sqrtq = sqrt((chi_ + x*invSig*x')*psi_);
    
    log_kernel = ...
        log(besselk( lambda_ - p/2, sqrtq )) ...
        - (p/2 - lambda_) * log( sqrtq );
    
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
    
    scoreSig = NaN(N,p_);
    score_nu = NaN(N,1);
    
    for ii = 1:N
        
        x = data_mat(ii,:) - mu_';
        sqrtq = sqrt((chi_ + x*invSig*x')*psi_);
       
        S = - invSigma ...
            - (besselk( lambda_ - p/2 - 1, sqrtq ) ...
                   + besselk( lambda_ - p/2 + 1, sqrtq )) ...
               / besselk( lambda_ - p/2, sqrtq ) / sqrtq ...
               * invSigma*x'*x*invSigma ...
            - (p/2 - lambda_)*invSigma*x'*x*invSigma/(chi_ + x*invSig*x');
        
        % Apply vec operator (and multiply with .5 factor).
        S = .5*S(:);
        
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
