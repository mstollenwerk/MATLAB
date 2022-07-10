function [ nLogL, logLcontr, varargout ] = matvitWishlike( Omega_, n, nu, X, varargin )
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
    Omega_ = ivechchol(all_param(1:p_));
    n = all_param(p_ + 1);
    nu = all_param(p_ + 2);
end
param.Omega_ = Omega_;
param.n = n;
param.nu = nu;
% This serves as a p.d. check on Omega_. Optionally this param struct can be returned.
param.all = [vechchol(Omega_); n; nu];
%% Log-Likelihood
logLcontr = NaN(N,1);

log_norm_const = ...
	gammaln( (n + nu*p)/2 ) ...
  - gammaln(n/2) ...
  - mvgammaln(nu/2, p) ...
  + n/2*log(n) ...
  + nu/2*log(det(Omega_));
	  
for ii = 1:N
    B = inv(X(:,:,ii));
 
    log_kernel = ...
		(nu+p+1)/2*logdet(B) ...
      - (n + nu*p)/2*log(n+trace(Omega_*B));
    
    logLcontr(ii) = log_norm_const + log_kernel;
end
nLogL = -sum(logLcontr);
%% Scores
if nargout >= 3
    invSig = inv(Omega_);
    
    score.Omega_ = NaN(N,p_);
    score.n = NaN(N,1);
    score.nu = NaN(N,1);
    
    for ii = 1:N
        
        B = inv(X(:,:,ii));
       
        trOmInvR = trace(Omega_*B);
        % General matrix derivative (ignoring symmetry of Omega_): [ enter to matrixcalculus.org: -dfn/2*log(det(Sigma)) - (dft + dfn*p)/2*log(dft+tr(inv(Sigma)*A)) ]        
        S = nu/2*invSig ...
            - (n + p*nu)/ 2 / (n + trOmInvR) * B;
        
        score.Omega_WishFishScaling(:,:,ii) = 2/nu*Omega_*S*Omega_;
        
        % Accounting for symmetry of Omega_:
        S = 2*S - diag(diag(S));
        
        % Numerical inaccuracies make S not exactly symmetric
        S = (S+S')./2;
                
        score.Omega_(ii,:) = vech(S);
        
        s_n = psi((n+p*nu)/2) -  nu*p/n - psi(n/2) - log(1+trOmInvR/n) + (n+p*nu)*(trOmInvR/n^2)/(1+trOmInvR/n);
        score.n(ii) = .5*s_n;
        score.n_scaled(ii) = -score.n(ii) ./ ...
            ( .25*psi(1,(n+p*nu)/2) - .25*psi(1,n/2) + .5*(p*nu+n+4)/(p*nu+n+2)*p*nu/(p*nu+n)/n );
        
        s_nu = p*psi((n+p*nu)/2) - p*log(n) - sum(mvpsi(ones(p,1)*nu/2)) - p*log(1+trOmInvR/n) + logdet(Omega_) - logdet(X(:,:,ii));
        score.nu(ii) = .5*s_nu;
        score.nu_scaled(ii) = -score.nu(ii) ./ ...
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
    fisherInfo.Omega_ = G'*(c1*kron2(invSig) + c2*vec2(invSig))*G;
    
    varargout{4} = fisherInfo;
    
end
%% Optional Parameter Vector Output
if nargout >= 5   
    varargout{3} = param;
end
end
