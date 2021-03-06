function [ nLogL, logLcontr, varargout ] = matvsWishlike( Sigma_, n, X, varargin )
%MATVSWISHLIKE Log-likelihood of Wishart distribution.
%
% Michael Stollenwerk
% michael.stollenwerk@live.com
% 30.08.2021

narginchk(3,4);
nargoutchk(0,6);

[p,~,N] = size(X);
p_ = p*(p+1)/2;
%% Parameters
if nargin == 4
    all_param = varargin{1};
    Sigma_ = ivechchol(all_param(1:p_));
    n = all_param(p_ + 1);
end

% Checking if Sigma_ is symmetric p.d.
param.Sigma_ = Sigma_;
param.n = n;
param.all = [vechchol(Sigma_); n];
%% Log-Likelihood
logLcontr = NaN(N,1);

log_norm_const = ...
    n*p/2*log(n) ...
  - n*p/2*log(2) ...
  - mvgammaln(n/2, p) ...
  - n/2*logdet(Sigma_);
	  
for ii = 1:N
    R = X(:,:,ii);
 
    log_kernel = ...
      -  trace(n*(Sigma_\R)) ...
      + (n-p-1)*logdet(R);
		
    logLcontr(ii) = log_norm_const + .5*log_kernel;
end
nLogL = -sum(logLcontr);
%% Scores
if nargout >= 3
    invSig = inv(Sigma_);    
    
    score.Sigma_ = NaN(p,p,N);
    
    for ii = 1:N
        
        R = X(:,:,ii);
        
        % General matrix derivative (ignoring symmetry of Sigma_): [ enter to matrixcalculus.org: -df/2*log(det(Sigma))-tr(inv(Sigma)*A) ]
        S = Sigma_\R/Sigma_ - invSig;
        S = n/2*S;
        
        score.SigmaNonSym = S;
        
        % Accounting for symmetry of Sigma_:
        S = S+S' - diag(diag(S));
        
        score.Sigma_(:,:,ii) = S;

        score.n(ii) = 1/2*(p+p*log(n/2)-sum(mvpsi(ones(p,1)*n/2))+logdet(R)-logdet(Sigma_)-trace(Sigma_\R));
        score.n_originalpdf(ii) = 1/2*(logdet(R)-p*log(2)-logdet(Sigma_/n)-sum(mvpsi(ones(p,1)*n/2)));
        score.n_originalpdf_scaled(ii) = 4*score.n_originalpdf(ii)/sum(mvpsi(ones(p,1)*n/2,1));
    
    end
    varargout{1} = score;
    
end
%% Hessian (Optional Output)

%% Optional Parameter Vector Output
if nargout >= 5   
    varargout{3} = param;
end
%% Fisher Info (Optional Output)
if nargout >= 6
    
    G = Dmatrix(p);
    
    fisherinfo.Sigma_ = n/2*G'*kron(invSig,invSig)*G;
    
    varargout{4} = fisherinfo;
end
end