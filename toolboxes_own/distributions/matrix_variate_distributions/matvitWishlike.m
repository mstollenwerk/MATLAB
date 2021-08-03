function [ nLogL, logLcontr, varargout ] = matvitWishlike( Sigma_, n, nu, X, varargin )
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
  + nu/2*log(nu) ...
  + n/2*log(det(Sigma_));
	  
for ii = 1:N
    B = inv(X(:,:,ii));
 
    log_kernel = ...
		(n+p+1)/2*logdet(B) ...
      - (nu + n*p)/2*log(nu+trace(Sigma_*B));
    
    logLcontr(ii) = log_norm_const + log_kernel;
end
nLogL = -sum(logLcontr);
%% Scores
if nargout >= 3
    invSig = inv(Sigma_);
    
    score.Sigma_ = NaN(N,p_);
    score.n = NaN(N,1);
    score.nu = NaN(N,1);
    
    for ii = 1:N
        
        B = inv(X(:,:,ii));
       
        % General matrix derivative (ignoring symmetry of Sigma_): [ enter to matrixcalculus.org: -dfn/2*log(det(Sigma)) - (dft + dfn*p)/2*log(dft+tr(inv(Sigma)*A)) ]        
        S = n/2*invSig ...
            - (nu + p*n)/ 2 / (nu + trace(Sigma_*B)) * B;
        
        % Accounting for symmetry of Sigma_:
        S = 2*S - diag(diag(S));
        
        % Numerical inaccuracies make S not exactly symmetric
        S = (S+S')./2;
                
        score.Sigma_(ii,:) = vech(S);

        % The score below are also easy to get quering wolframalpha.com 
        % with eg "d/da (log(Gamma(1/2 (a+ p n))))".
        % I am just too lazy to write them down right now.
        score.n(ii) = NaN;
        score.nu(ii) = NaN;
    
    end
    
    varargout{1} = score;
    
end
%% Hessian (Optional Output)

%% Fisher Info (Optional Output)
% enter to matrixcalculus.org: -(dfn/2*inv(Sigma)-(dft+dfn*p)/(2*(dft+tr(A*inv(Sigma))))*inv(Sigma)*A*inv(Sigma)), change the sign, then take expectation to arrive at
if nargout >= 6
    
    G = Dmatrix(p);
    c1 = n_/2*(nu_+p*n_)/(nu_+p*n_+2);
    c2 = -n_^2/2/(nu_+p*n_+2);
    fisherInfo.Sigma_ = -G'*(c1*kron2(invSig) + c2*vec2(invSig))*G;
    
    varargout{2} = fisherInfo;
    
end
%% Optional Parameter Vector Output
if nargout >= 5   
    varargout{3} = param;
end
end
