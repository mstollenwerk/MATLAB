function [ nLogL, logLcontr, varargout ] = matviWishlike( Omega_, nu, R, varargin )
%MATVIWISHLIKE Log-likelihood of inverse-Wishart matrix.
%
% Michael Stollenwerk
% michael.stollenwerk@live.com

narginchk(3,4);
nargoutchk(0,6);

[p,~,N] = size(R);
p_ = p*(p+1)/2;
%% Parameters
if nargin == 4
    all_param = varargin{1};
    Omega_ = ivechchol(all_param(1:p_));
    nu = all_param(p_ + 1);
end
% Checking if Omega_ is symmetric p.d.
param.Omega_ = Omega_;
param.n = nu;
param.all = [vechchol(Omega_); nu];
%% Log-Likelihood
logLcontr = NaN(N,1);

log_norm_const = ...
  - nu*p/2*log(2) ...
  - mvgammaln(nu/2, p) ...
  + nu/2*logdet(Omega_);
	  
for ii = 1:N
    B = inv(R(:,:,ii));
 
    log_kernel = ...
      - trace(Omega_*B) ...
      + (nu+p+1)*logdet(B);
		
    logLcontr(ii) = log_norm_const + .5*log_kernel;
end
nLogL = -sum(logLcontr);
%% Scores
if nargout >= 3
    invSig = inv(Omega_);    
    
    score.Omega_ = NaN(N,p_);
    score.Omega_WishFishScaling = NaN(p,p,N);
    score.nu = NaN(N,1);
    for ii = 1:N
        
        invA = inv(R(:,:,ii));
       
        % General matrix derivative (ignoring symmetry of Omega_):
        S = .5*( nu*invSig - invA );
        
        score.Omega_WishFishScaling(:,:,ii) = 2/nu*Omega_*S*Omega_;

        % Accounting for symmetry of Omega_:
        S = 2*S - diag(diag(S));
        
        % Numerical inaccuracies make S not exactly symmetric
        S = (S+S')./2;
        
        score.Omega_(ii,:) = vech(S);

        score.nu(ii) = .5*( -p*log(2) - sum(mvpsi(ones(p,1)*nu/2)) + logdet(Omega_) - logdet(R(:,:,ii)) );
        
        score.nu_scaled(ii) = - score.nu(ii) ./ ...
            (-.25*sum(mvpsi(ones(p,1)*nu/2,1)));
    
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
    fisherinfo.Omega_ = nu/2*G'*kron(invSig,invSig)*G;
    
    fisherinfo.df = NaN;
    
    varargout{4} = fisherinfo;
    
end
end
