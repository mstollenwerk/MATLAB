function [ nLogL, logLcontr, varargout ] = matvWishlike( Omega_, n, R, varargin )
%MATVWISHLIKE Log-likelihood of Wishart distribution.
%
% Michael Stollenwerk
% michael.stollenwerk@live.com
% 14.08.2020

narginchk(3,4);
nargoutchk(0,6);

[p,~,T] = size(R);
p_ = p*(p+1)/2;
%% Parameters
if nargin == 4
    all_param = varargin{1};
    Omega_ = ivechchol(all_param(1:p_));
    n = all_param(p_ + 1);
end

% Checking if Omega_ is symmetric p.d.
param.Omega_ = Omega_;
param.n = n;
param.all = [vechchol(Omega_); n];
%% Log-Likelihood
logLcontr = NaN(T,1);

log_norm_const = ...
  - n*p/2*log(2) ...
  - mvgammaln(n/2, p) ...
  - n/2*logdet(Omega_);
	  
for tt = 1:T
    A = R(:,:,tt);
 
    log_kernel = ...
      -  trace(Omega_\A) ...
      + (n-p-1)*logdet(A);
		
    logLcontr(tt) = log_norm_const + .5*log_kernel;
end
nLogL = -sum(logLcontr);
%% Scores
if nargout >= 3
    invSig = inv(Omega_);    
    
    score.Omega_ = NaN(T,p_);
    score.Omega_WishFishScaling = NaN(p,p,T);
    score.df = NaN(T,1);
    
    for tt = 1:T
        
        A = R(:,:,tt);
        
        % General matrix derivative (ignoring symmetry of Omega_): [ enter to matrixcalculus.org: -df/2*log(det(Sigma))-tr(inv(Sigma)*A) ]
        S = - n/2*invSig ...
            + 1/2*(Omega_\A/Omega_);
        
        score.Omega_WishFishScaling(:,:,tt) = 2/n*Omega_*S*Omega_;
        
        % Accounting for symmetry of Omega_:
        S = 2*S - diag(diag(S));
        
        % Numerical inaccuracies make S not exactly symmetric
        S = (S+S')./2;
        
        score.Omega_(tt,:) = vech(S);

        score.n(tt) = 1/2*(logdet(A)-p*log(2)-logdet(Omega_)-sum(mvpsi(ones(p,1)*n/2)));
        score.n_scaled(tt) = 4*score.n(tt)/sum(mvpsi(ones(p,1)*n/2,1));
    
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
    
    fisherinfo.Omega_ = n/2*G'*kron(invSig,invSig)*G;
    
    varargout{4} = fisherinfo;
end
end
