function [ nLogL, logLcontr, varargout ] = ...
    matvitRiesz2like( Omega_, n, nu, X, varargin )
%MATVtRieszLIKE Negative log-likelihood and score of
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
if nargin == 5 %%%%%%%
    if ~(isempty(Omega_) && isempty(n) && isempty(nu))
        error('Cannot input all_param and any of the parameters individually!')
    end
    all_param = varargin{1};
    Omega_ = ivechchol(all_param(1 : p_));
    n = all_param(p_ +  1 : p_ + p);
    nu = all_param(p_ + p + 1);    
end
% Checking if Omega_ is symmetric p.d.
param.Omega_ = Omega_;
C = chol(Omega_,'lower');
param.chol_Omega_ = vech(C);
param.nu = nu;
param.n = n;
param.all = [param.chol_Omega_; n; nu];
%% Log-likelihood computation
logLcontr = NaN(N,1);

term1 = gammaln((n + sum(nu))/2);
term2 = -gammaln(n/2);
term3 = -ugmvgammaln(nu./2);
term4 = -sum(nu)/2*log(n);
term5 = loglpwdet([],nu./2,diag(C)); % upwdet(invS,-n) = lpwdet(S,n)

log_normalizing_constant = term1 + term2 + term3 + term4 + term5;

for ii = 1:N
    
    A = X(:,:,ii);
    
    term6 = loglpwdet(A,-(nu-p-1)./2); % upwdet(invS,-n) = lpwdet(S,n)
    term7 = -(n + sum(nu))/2*log(1 + trace(Omega_/A)/n);
    term8 = -(p+1)*logdet(A);
    
    log_kernel = term6 + term7 + term8;

    logLcontr(ii) = log_normalizing_constant + log_kernel;
end
nLogL = -sum(logLcontr);
%% Score computation
if nargout >= 3
    
    score.Omega_ = NaN(N,p_);
    score.n = NaN(N,p);
    score.nu = NaN(N,1);
    
    q = C'\diag(nu)/C;
    
    for ii = 1:N
        
        A = X(:,:,ii);
        
        % General matrix derivative (ignoring symmetry of Omega_):
        S = q - (n + sum(nu))/(n + trace(Omega_/A))*inv(A);
        S = .5*S;
        
        score.Omega_WishFishScaling(:,:,ii) = 2/mean(nu)*Omega_*S*Omega_;
        
        % Accounting for symmetry of Omega_:
        S = 2*S - diag(diag(S));
        
        % Numerical inaccuracies make S not exactly symmetric
        S = (S+S')./2;
                
        score.Omega_(ii,:) = vech(S);
        
        score.n(ii,:) = NaN(p,1);
        score.nu(ii) = NaN;
    end
    
    varargout{1} = score;

end
%% Hessian (Optional Output)

%% Fisher Info (Optional Output)
if nargout >= 6
    
    [G,iG] = Dmatrix(p);
    I = speye(p);
    L = ELmatrix(p);
    invSig = inv(Omega_);
    
    term1 = G'*kron(inv(C)',q)*L'/(iG*kron(C,I)*L');
    term2 = -1/(n+sum(nu)+2)*G'*(2*kron(invSig,q) + vec2(q))*G;
    
    fisherinfo.Omega_ = .5*(term1+term2);
    fisherinfo.n = NaN;
    
    varargout{4} = fisherinfo;
end
%% Optional Parameter Vector Output
if nargout >= 5
    varargout{3} = param;
end
end