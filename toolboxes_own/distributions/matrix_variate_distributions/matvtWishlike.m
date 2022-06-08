function [ nLogL, logLcontr, varargout ] = matvtWishlike( Omega_, n, nu, R, varargin )
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
%   Omega_.
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
    if ~(isempty(Omega_) && isempty(n) && isempty(nu))
        error('Cannot input all_param and any of the parameters individually!')
    end
    all_param = varargin{1};
    Omega_ = ivechchol(all_param(1:p_));
    n = all_param(p_ + 1);
    nu = all_param(p_ + 2);
end
% Checking if Sigma_ is symmetric p.d.
param.Omega_ = Omega_;
param.chol_Omega_ = vechchol(Omega_);
param.n = n;
param.nu = nu;
param.all = [param.chol_Omega_; n; nu];
%% Log-Likelihood
logLcontr = NaN(N,1);
trQ = NaN(N,1);

for ii = 1:N
    A = R(:,:,ii);
    trQ(ii) = trace(Omega_\A);

    term1 = gammaln( (nu + n*p)/2 );
    term2 = - gammaln(nu/2);
    term3 = - mvgammaln(n/2, p);
    term4 = nu/2*log(nu);
    term5 = - n/2*logdet(Omega_);    
    log_norm_const = term1 + term2 + term3 + term4 + term5;
    
    term6 = (n-p-1)/2*logdet(A);
    term7 = - (nu + n*p)/2*log(nu+trQ(ii));
    log_kernel = term6 + term7;    
    
    logLcontr(ii) = log_norm_const + log_kernel;
end
nLogL = -sum(logLcontr);
%% Scores (Optional Output)
if nargout >= 3
    
    score.Omega_ = NaN(N,p_);
    score.Omega_WishFishScaling = NaN(p,p,N);
    score.cholSigma = NaN;
    score.n = NaN;
    score.nu = NaN;
    score.all = NaN;
    
    invSig = inv(Omega_);    
    for ii = 1:N        
        A = R(:,:,ii);
       
        % General matrix derivative (ignoring symmetry of Omega_): [ enter to matrixcalculus.org: -dfn/2*log(det(Sigma)) - (dft + dfn*p)/2*log(dft+tr(inv(Sigma)*A)) ]        
        S = - n*invSig ...
            + (nu + p*n) / (nu + trace(Omega_\A)) * (Omega_\A/Omega_);
        S = .5*S;
        
        score.Omega_WishFishScaling(:,:,ii) = 2/n*Omega_*S*Omega_;
        
        % Accounting for symmetry of Omega_:
        S = 2*S - diag(diag(S));
        
        % Numerical inaccuracies make S not exactly symmetric
        S = vech((S+S')./2);

        score.Omega_(ii,:) = S;

        % Wolframalpha querie: d/da log(gamma((a + p*n)/2)) - p*n/2*log(a)
        % - log(gamma(a/2)) - (a + p*n)/2*log(1+q/a) 
        term1 = trQ(ii)*(nu + n*p)/nu^2/(trQ(ii)/nu+1);
        term2 = -n*p/nu;
        term3 = psi(.5*(nu+n*p));
        term4 = log(trQ(ii)/nu+1);
        term5 = psi(nu/2);
        score.nu(ii) = .5*(term1+term2+term3+term4+term5);
        
        score.n(ii) = NaN;
    
    end
    
    varargout{1} = score;
    
end
%% Hessian (Optional Output)
if nargout >= 4
    
    hessian.Omega_ = NaN;
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
    fisherInfo.Omega_ = G'*(c1*kron2(invSig) + c2*vec2(invSig))*G;
    
    varargout{4} = fisherInfo;
end
%% Optional Parameter Vector Output
if nargout >= 5    
    varargout{3} = param;
end
end
