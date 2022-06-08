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
    score.nu = NaN;
    
    invSig = inv(Sigma_);    
    for ii = 1:N     
        
        R_ = R(:,:,ii);
       
        % General matrix derivative (ignoring symmetry of Sigma_): [ enter to matrixcalculus.org: -dfn/2*log(det(Sigma)) - (dft + dfn*p)/2*log(dft+tr(inv(Sigma)*A)) ]        
        S = (nu+p*n)/(nu-2+n*trace(Sigma_\R_)) * (Sigma_\R_/Sigma_) - invSig;
            
        S = n/2*S;
        
        score.SigmaNonSym = S;
        
        % Accounting for symmetry of Sigma_:
        S = S+S' - diag(diag(S));

        score.Sigma_(:,:,ii) = S;

        % Wolframalpha querie: d/da log(gamma((a + p*n)/2)) - p*n/2*log(a)
        % - log(gamma(a/2)) - (a + p*n)/2*log(1+q/a) 
        term1 = trQ(ii)*(nu + n*p)/nu^2/(trQ(ii)/nu+1);
        term2 = -n*p/nu;
        term3 = psi(.5*(nu+n*p));
        term4 = log(trQ(ii)/nu+1);
        term5 = psi(nu/2);
        score.nu(ii) = .5*(term1+term2+term3+term4+term5);
    
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
