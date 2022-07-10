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
    nu = all_param(p_ + 1);
    n = all_param(p_ + 2);
end
param.Sigma_ = Sigma_;
param.nu = nu;
param.n = n;
% This serves as a p.d. check on Sigma_. Optionally this param struct can be returned.
param.all = [vechchol(Sigma_); n; nu];
%% Log-Likelihood
logLcontr = NaN(N,1);

log_norm_const = ...
	gammaln( (n + nu*p)/2 ) ...
  - gammaln(n/2) ...
  - mvgammaln(nu/2, p) ...
  + nu*p/2*log((nu-p-1)/n) ...
  + nu/2*logdet(Sigma_);
	  
for ii = 1:N
    R = X(:,:,ii);
 
    log_kernel = ...
      - (nu+p+1)/2*logdet(R) ...
      - (n + nu*p)/2*log(1+(nu-p-1)/n*trace(Sigma_/R));
    
    logLcontr(ii) = log_norm_const + log_kernel;
end
nLogL = -sum(logLcontr);
%% Scores
if nargout >= 3
    invSig = inv(Sigma_);
    
    score.Sigma_ = NaN(p,p,N);
    
    for ii = 1:N
        
        R = X(:,:,ii);
        
        trSigInvR = trace(Sigma_/R);
        % General matrix derivative (ignoring symmetry of Sigma_): [ enter to matrixcalculus.org: -dfn/2*log(det(Sigma)) - (dft + dfn*p)/2*log(dft+tr(inv(Sigma)*A)) ]        
        S = nu*invSig - (n+p*nu)/(n/(nu-p-1) + trSigInvR)*inv(R);
        S = .5.*S;
        
        score.SigmaNonSym = S;
        
        % Accounting for symmetry of Sigma_:
        S = S+S' - diag(diag(S));
                
        score.Sigma_(:,:,ii) = S;
        
        trOmInvR = (nu-p-1)*trSigInvR;
        s_n = psi((n+p*nu)/2) -  nu*p/n - psi(n/2) - log(1+trOmInvR/n) + (n+p*nu)*(trOmInvR/n^2)/(1+trOmInvR/n);
        score.n_originalpdf(ii) = .5*s_n;
        
        s_nu = p*psi((n+p*nu)/2) - p*log(n) - sum(mvpsi(ones(p,1)*nu/2)) - p*log(1+trOmInvR/n) + logdet((nu-p-1)*Sigma_) - logdet(R);
        score.nu_originalpdf(ii) = .5*s_nu;
        
        score.n_originalpdf_scaled(ii) = -score.n_originalpdf(ii) ./ ...
            ( .25*psi(1,(n+p*nu)/2) - .25*psi(1,n/2) + .5*(p*nu+n+4)/(p*nu+n+2)*p*nu/(p*nu+n)/n );

        score.nu_originalpdf_scaled(ii) = -score.nu_originalpdf(ii) ./ ...
            ( .25*( p^2*psi(1,(n+p*nu)/2) - sum(mvpsi(ones(p,1)*nu/2,1)) ) );        
        
    end
    
    varargout{1} = score;
    
end
%% Hessian (Optional Output)

%% Fisher Info (Optional Output)
% enter to matrixcalculus.org: -(dfn/2*inv(Sigma)-(dft+dfn*p)/(2*(dft+tr(A*inv(Sigma))))*inv(Sigma)*A*inv(Sigma)), change the sign, then take expectation to arrive at
if nargout >= 6
    
    G = Dmatrix(p);
    c1 = nu/2*(n+p*nu)/(n+p*nu+2);
    c2 = -nu^2/2/(n+p*nu+2);
    fisherInfo.Sigma_ = G'*(c1*kron2(invSig) + c2*vec2(invSig))*G;
    
    varargout{4} = fisherInfo;
    
end
%% Optional Parameter Vector Output
if nargout >= 5   
    varargout{3} = param;
end
end
