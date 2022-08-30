function [ nLogL, logLcontr, varargout ] = matvsiWishlike( Sigma_, nu, R, varargin )
%MATVSIWISHLIKE Log-likelihood of standardized inverse-Wishart matrix.
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
    Sigma_ = ivechchol(all_param(1:p_));
    nu = all_param(p_ + 1);
end
% Checking if Sigma_ is symmetric p.d.
param.Sigma_ = Sigma_;
param.nu = nu;
param.all = [vechchol(Sigma_); nu];
%% Log-Likelihood
logLcontr = NaN(N,1);

log_norm_const = ...
    nu*p/2*log((nu-p-1)/2) ...
  - mvgammaln(nu/2, p) ...
  + nu/2*logdet(Sigma_);
	  
for ii = 1:N
    
    R_ = R(:,:,ii);
 
    log_kernel = ...
      - (nu-p-1)*trace(Sigma_/R_) ...
      - (nu+p+1)*logdet(R_);
		
    logLcontr(ii) = log_norm_const + .5*log_kernel;
end
nLogL = -sum(logLcontr);
%% Scores
if nargout >= 3
    invSig = inv(Sigma_);    
    
    score.Sigma_ = NaN(p,p,N);
    score.nu_originalpdf = NaN(N,1);
    
    for ii = 1:N
        
        invR = inv(R(:,:,ii));
       
        % General matrix derivative (ignoring symmetry of Sigma_):
        S = .5*( nu*invSig - (nu-p-1)*invR );
        
        score.SigmaNonSym = S;
        
        % Accounting for symmetry of Sigma_:
        S = S+S' - diag(diag(S));
        
        score.Sigma_(:,:,ii) = S;
            
        Omega_ = (nu-p-1)*Sigma_;
        score.nu_originalpdf(ii) = ...
            .5*( -p*log(2) - sum(mvpsi(ones(p,1)*nu/2)) + logdet(Omega_) - logdet(R(:,:,ii)) );
    
        score.nu_originalpdf_scaled(ii) = - score.nu_originalpdf(ii) ./ ...
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
    fisherinfo.Sigma_ = nu/2*G'*kron(invSig,invSig)*G;
    
    fisherinfo.df = NaN;
    
    varargout{4} = fisherinfo;
    
end
end
