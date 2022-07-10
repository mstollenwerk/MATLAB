function [ nLogL, logLcontr, varargout ] = matvstWishlike( Sigma_, n, nu, R, varargin )
%MATVTWISHLIKE
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
%   This is the distribution of the quadratic form XX', where the columns
%   of the p by df_n matrix X follow multivariate central t distributions
%   with the same degree of freedom df_t and the same dispersion matrix 
%   Sigma_.
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

[p,~,N] = size(R);
p_ = p*(p+1)/2;
%% Parameters
if nargin >= 5
    if ~(isempty(Sigma_) && isempty(n) && isempty(nu))
        error('Cannot input all_param and any of the parameters individually!')
    end
    all_param = varargin{1};
    Sigma_ = ivechchol(all_param(1:p_));
    n = all_param(p_ + 1);
    nu = all_param(p_ + 2);
end
% Checking if Sigma_ is symmetric p.d.
param.Sigma_ = Sigma_;
param.chol_Sigma_ = vechchol(Sigma_);
param.n = n;
param.nu = nu;
param.all = [param.chol_Sigma_; n; nu];
%% Log-Likelihood
logLcontr = NaN(N,1);
trQ = NaN(N,1);

term1 = gammaln( (nu + n*p)/2 );
term2 = - gammaln(nu/2);
term3 = - mvgammaln(n/2, p);
term4 = n*p/2*log(n/(nu-2));
term5 = - n/2*logdet(Sigma_);

log_norm_const = term1 + term2 + term3 + term4 + term5;

for ii = 1:N
    
    R_ = R(:,:,ii);
    trQ(ii) = trace(Sigma_\R_);
    
    term6 = (n-p-1)/2*logdet(R_);
    term7 = - (nu + n*p)/2*log(1+n/(nu-2)*trQ(ii));
    log_kernel = term6 + term7;    
    
    logLcontr(ii) = log_norm_const + log_kernel;
end
nLogL = -sum(logLcontr);
%% Scores (Optional Output)
if nargout >= 3
    
    score.Sigma_ = NaN(p,p,N);
    score.nu_originalpdf = NaN(N,1);
    score.n_originalpdf = NaN(N,1);
    
    invSig = inv(Sigma_);    
    for ii = 1:N     
        
        R_ = R(:,:,ii);
       
        trInvSigR = trace(Sigma_\R_);
        % General matrix derivative (ignoring symmetry of Sigma_): [ enter to matrixcalculus.org: -dfn/2*log(det(Sigma)) - (dft + dfn*p)/2*log(dft+tr(inv(Sigma)*A)) ]        
        S = (nu+p*n)/(nu-2+n*trInvSigR) * (Sigma_\R_/Sigma_) - invSig;
            
        S = n/2*S;
        
        score.SigmaNonSym = S;
        
        % Accounting for symmetry of Sigma_:
        S = S+S' - diag(diag(S));

        score.Sigma_(:,:,ii) = S;
        
        score.n(ii) = .5*( p*log(n/(nu-2)) + p + p*psi((nu+p*n)/2) - sum(mvpsi(ones(p,1)*n/2)) ...
                          - logdet(Sigma_) + logdet(R_) - p*log(1+trInvSigR*n/(nu-2)) ...
                          - (nu+p*n)/2*trInvSigR/(nu-2)/(1+trInvSigR*n/(nu-2)) );
                      
        score.nu(ii) = .5*( -p*n/(nu-2) + psi((nu+p*n)/2) - psi(nu/2) - p*n/nu ...
                         - log(1+trInvSigR*n/(nu-2)) + (nu+p*n)*n/(nu-2)^2*trInvSigR/(1+trInvSigR*n/(nu-2)) );                      
        
        Omega_ = Sigma_*(nu-2)/nu/n;
        trInvOmR = trInvSigR/(nu-2)*nu*n;
        score.n_originalpdf(ii) = .5*( p*psi((nu+p*n)/2) - sum(mvpsi(ones(p,1)*n/2)) ...
                          - p*log(nu) - logdet(Omega_) + logdet(R_) - p*log(1+trInvOmR/nu) );
        
        score.nu_originalpdf(ii) = .5*( psi((nu+p*n)/2) - psi(nu/2) - p*n/nu ...
                         - log(1+trInvOmR/nu) + (nu+p*n)*(trInvOmR/nu^2)/(1+trInvOmR/nu) );
                     
        score.n_originalpdf_scaled(ii) = -score.n_originalpdf(ii) / ...
            (.25*p^2*psi(1,(nu+p*n)/2) - .25*sum(mvpsi(ones(p,1)*n/2,1)));
                   
        score.nu_originalpdf_scaled(ii) = -score.nu_originalpdf(ii) / ...
            (.25*psi(1,(nu+p*n)/2) - .25*psi(1,nu/2) + .5*(p*n+nu+4)/(p*n+nu+2)*p*n/(p*n+nu)/nu);                        
    
    end
    
    varargout{1} = score;
    
end
%% Hessian (Optional Output)
if nargout >= 4
    
    hessian.Sigma_ = NaN;
    hessian.cholSigma = NaN;
    hessian.n = NaN;
    hessian.nu = NaN;
    hessian.all = NaN;    
    
    varargout{2} = hessian;
end
%% Fisher Info (Optional Output)
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
