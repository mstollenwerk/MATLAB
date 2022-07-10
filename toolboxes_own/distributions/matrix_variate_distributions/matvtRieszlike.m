function [ nLogL, logLcontr, varargout ] = ...
    matvtRieszlike( Omega_, n, nu, X, varargin )
%MATVRIESZTLIKE Negative log-likelihood and score of
%
% USAGE:
%   [ nLogL, score, param ] = matvFlike( Omega_, df_1, df_2, X )
%   [ nLogL, score, param ] = matvFlike( [], [], [], X, all_param )
%
% INPUTS:
%   X       - Array (p by p). Symmetric p.d. data matrix.
%   Omega_  - Array (p by p). Symmetric p.d. parameter matrix, 
%             regulates covariance. 
%   DF_1    - Double. First degrees of freedom parameter. 
%   DF_2    - Double. Second degrees of freedom parameter. 
%
% OUTPUTS:
%   NLOGL   - Double. Negative log-likelihood value.
%   SCORE   - Struct. Fields as parameter names. Contains derivatives
%             of log-likelihood w.r.t. parameters.
%
% See also MVBETALN MATVFRND IVECHCHOL
% 
% COMMENTS:
%
%   Derivatives were derived using matrixcalculus.org
%   References [3] and [4] use slightly different definitions of the
%   distributions.
%
% REFERENCES:
%   [1] Gupta and Nagar (2001) - Matrix Variate Distributions, p.156.
%   [2] Mulder and Pericchi (2018) - The Matrix F Prior for Estimating and
%           Testing Covariance Matrices.
%   [3] Opschoor, Janus, Lucas and Van Dijk (2018) - New HEAVY Models for 
%           Fat-Tailed Realized Covariances and Returns.
%   [4] Bodnar, Okhrin, Parolya and Stollenwerk (2020)
%
% DEPENDENCIES:
%   MVBETALN IVECHCHOL
%
% Michael Stollenwerk
% michael.stollenwerk@live.com
% 11.02.2021

[p,~,N] = size(X);
p_ = p*(p+1)/2;
narginchk(4,5); %%%%%%%
nargoutchk(0,6);
%% Param
if nargin >= 5 %%%%%%%
    if ~(isempty(Omega_) && isempty(nu) && isempty(n))
        error('Cannot input all_param and any of the parameters individually!')
    end
    all_param = varargin{1};
    Omega_ = ivechchol(all_param(1 : p_));
    n = all_param(p_ + 1 : p_ + p);
    nu = all_param(p_ + p + 1);
end
% Checking if Omega_ is symmetric p.d.
param.Omega_ = Omega_;
COm = chol(Omega_,'lower');
param.chol_Omega_ = vech(COm);
param.nu = nu;
param.n = n;
param.all = [param.chol_Omega_; n; nu];
%% Log-likelihood computation
logLcontr = NaN(N,1);

term1 = gammaln((nu + sum(n))/2);
term2 = -gammaln(nu/2);
term3 = -lgmvgammaln(n./2);
term4 = -sum(n)/2*log(nu);
term5 = -loglpwdet([],n./2,diag(COm));

log_normalizing_constant = term1 + term2 + term3 + term4 + term5;

for ii = 1:N
    term6 = loglpwdet(X(:,:,ii),(n-p-1)./2);
    term7 = -(nu + sum(n))/2*log(1 + trace(Omega_\X(:,:,ii))/nu);
    
    log_kernel = term6 + term7;

    logLcontr(ii) = log_normalizing_constant + log_kernel;
end
nLogL = -sum(logLcontr);
%% Score computation
if nargout >= 3
    
    score.Omega_ = NaN(N,p_);
    score.n = NaN(N,p);
    score.nu = NaN(N,1);
    
    for ii = 1:N
        
        A = X(:,:,ii);
        
        % General matrix derivative (ignoring symmetry of Omega_):
        invOmA = Omega_\A;
        S = (nu+sum(n))/(nu+trace(invOmA))*invOmA/Omega_ - COm'\diag(n)/COm;
        S = .5*S;
        
        score.Omega_WishFishScaling(:,:,ii) = 2/mean(n)*Omega_*S*Omega_;
        
        % Accounting for symmetry of Omega_:
        S = 2*S - diag(diag(S));
        
        % Numerical inaccuracies make S not exactly symmetric
        S = (S+S')./2;
                
        score.Omega_(ii,:) = vech(S);
        
        trInvOmA = trace(invOmA);
        score.n(ii,:) = ...
            1/2*( psi((nu+sum(n))/2) - mvpsi(n/2) - log(nu + trInvOmA) ) ...
          - log(diag(COm)) + log(diag(chol(A,'lower')));
        score.nu(ii) = ...
            1/2*( psi((nu+sum(n))/2) - psi(nu/2) ...
          - log(1 + trInvOmA/nu) - (sum(n)-trInvOmA)/(nu+trInvOmA) );
        
        score.n_scaled(ii,:) = -score.n(ii,:)' ./ ...
            (.25*psi(1,(nu+sum(n))/2) - .25*mvpsi(n/2,1));
                   
        score.nu_scaled(ii) = -score.nu(ii) ./ ...
            (.25*psi(1,(nu+sum(n))/2) - .25*psi(1,nu/2) + .5*(sum(n)+nu+4)/(sum(n)+nu+2)*sum(n)/(sum(n)+nu)/nu);
        
    end
    
    varargout{1} = score;

end
%% Hessian (Optional Output)

%% Fisher Info (Optional Output)
if nargout >= 6
    
    G = Dmatrix(p);
    I = speye(p);
    L = ELmatrix(p);
    
    invSig = inv(Omega_);
    q = COm'\diag(n)/COm;
    
    term1 = kron(COm'\diag(n),invSig)*L'/(G'*kron(COm,I)*L')*G';
    term2 = -1/(nu+sum(n)+2)*(2*kron(q,invSig) + vec2(q));
    
    fisherinfo.Omega_ = .5*G'*(term1+term2)*G;
    fisherinfo.n = NaN;
    
    varargout{4} = fisherinfo;
end
%% Optional Parameter Vector Output
if nargout >= 5   
    varargout{3} = param;
end
end