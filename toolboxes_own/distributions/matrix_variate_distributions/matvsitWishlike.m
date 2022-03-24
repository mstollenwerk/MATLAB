function [ nLogL, logLcontr, varargout ] = matvsitWishlike( Sigma_, n, nu, X, varargin )
%MATVITWISHLIKE
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
% 11.08.2020
%
% DEPENDENCIES:

narginchk(4,5);
nargoutchk(0,6);

[p,~,N] = size(X);
p_ = p*(p+1)/2;
%% Parameters
if nargin >= 5
    all_param = varargin{1};
    Sigma_ = ivechchol(all_param(1:p_));
    n = all_param(p_ + 1);
    nu = all_param(p_ + 2);
end
param.Sigma_ = Sigma_;
param.n = n;
param.nu = nu;
% This serves as a p.d. check on Sigma_. Optionally this param struct can be returned.
param.all = [vechchol(Sigma_); n; nu];
%% Log-Likelihood
logLcontr = NaN(N,1);

log_norm_const = ...
	gammaln( (nu + n*p)/2 ) ...
  - gammaln(nu/2) ...
  - mvgammaln(n/2, p) ...
  + n*p/2*log((n-p-1)/nu) ...
  + n/2*logdet(Sigma_);
	  
for ii = 1:N
    R = X(:,:,ii);
 
    log_kernel = ...
      - (n+p+1)/2*logdet(R) ...
      - (nu + n*p)/2*log(1+(n-p-1)/nu*trace(Sigma_/R));
    
    logLcontr(ii) = log_norm_const + log_kernel;
end
nLogL = -sum(logLcontr);
%% Scores
if nargout >= 3
    invSig = inv(Sigma_);
    
    score.Sigma_ = NaN(p,p,N);
    
    for ii = 1:N
        
        R = X(:,:,ii);
       
        % General matrix derivative (ignoring symmetry of Sigma_): [ enter to matrixcalculus.org: -dfn/2*log(det(Sigma)) - (dft + dfn*p)/2*log(dft+tr(inv(Sigma)*A)) ]        
        S = n*invSig - (nu+p*n)/(nu/(n-p-1) + trace(Sigma_/R))*inv(R);
        S = .5.*S;
        
        % Accounting for symmetry of Sigma_:
        S = S+S' - diag(diag(S));
                
        score.Sigma_(:,:,ii) = S;
    
    end
    
    varargout{1} = score;
    
end
%% Hessian (Optional Output)

%% Fisher Info (Optional Output)
% enter to matrixcalculus.org: -(dfn/2*inv(Sigma)-(dft+dfn*p)/(2*(dft+tr(A*inv(Sigma))))*inv(Sigma)*A*inv(Sigma)), change the sign, then take expectation to arrive at
if nargout >= 6
    
    G = Dmatrix(p);
    c1 = n/2*(nu+p*n)/(nu+p*n+2);
    c2 = -n^2/2/(nu+p*n+2);
    fisherInfo.Sigma_ = G'*(c1*kron2(invSig) + c2*vec2(invSig))*G;
    
    varargout{4} = fisherInfo;
    
end
%% Optional Parameter Vector Output
if nargout >= 5   
    varargout{3} = param;
end
end
