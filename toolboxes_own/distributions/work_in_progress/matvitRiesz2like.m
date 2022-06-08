function [ nLogL, logLcontr, varargout ] = ...
    matvitRiesz2like( Sigma_, n, nu, X, varargin )
%MATVtRieszLIKE Negative log-likelihood and score of
%
% USAGE:
%   [ nLogL, score, param ] = matvFlike( Sigma_, df_1, df_2, X )
%   [ nLogL, score, param ] = matvFlike( [], [], [], X, all_param )
%
% INPUTS:
%   X       - Array (p by p). Symmetric p.d. data matrix.
%   SIGMA_  - Array (p by p). Symmetric p.d. parameter matrix, 
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
    if ~(isempty(Sigma_) && isempty(nu) && isempty(n))
        error('Cannot input all_param and any of the parameters individually!')
    end
    all_param = varargin{1};
    Sigma_ = ivechchol(all_param(1 : p_));
    n = all_param(p_ +  1 : p_ + p);
    nu = all_param(p_ + p + 1);    
end
% Checking if Sigma_ is symmetric p.d.
param.Sigma_ = Sigma_;
C = chol(Sigma_,'lower');
param.chol_Sigma_ = vech(C);
param.nu = nu;
param.n = n;
param.all = [param.chol_Sigma_; n; nu];
%% Log-likelihood computation
logLcontr = NaN(N,1);

term1 = gammaln((nu + sum(n))/2);
term2 = -gammaln(nu/2);
term3 = -ugmvgammaln(n./2);
term4 = -sum(n)/2*log(nu);
term5 = loglpwdet([],n./2,diag(C)); % upwdet(invS,-n) = lpwdet(S,n)

log_normalizing_constant = term1 + term2 + term3 + term4 + term5;

for ii = 1:N
    
    A = X(:,:,ii);
    
    term6 = loglpwdet(A,-(n-p-1)./2); % upwdet(invS,-n) = lpwdet(S,n)
    term7 = -(nu + sum(n))/2*log(1 + trace(Sigma_/A)/nu);
    term8 = -(p+1)*logdet(A);
    
    log_kernel = term6 + term7 + term8;

    logLcontr(ii) = log_normalizing_constant + log_kernel;
end
nLogL = -sum(logLcontr);
%% Score computation
if nargout >= 3
    
    score.Sigma_ = NaN(N,p_);
    score.n = NaN(N,p);
    score.nu = NaN(N,1);
    
    q = C'\diag(n)/C;
    
    for ii = 1:N
        
        A = X(:,:,ii);
        
        % General matrix derivative (ignoring symmetry of Sigma_):
        S = q - (nu + sum(n))/(nu + trace(Sigma_/A))*inv(A);
        S = .5*S;
        
        score.Sigma_WishFishScaling(:,:,ii) = 2/mean(n)*Sigma_*S*Sigma_;
        
        % Accounting for symmetry of Sigma_:
        S = 2*S - diag(diag(S));
        
        % Numerical inaccuracies make S not exactly symmetric
        S = (S+S')./2;
                
        score.Sigma_(ii,:) = vech(S);
        
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
    invSig = inv(Sigma_);
    
    term1 = G'*kron(inv(C)',q)*L'/(iG*kron(C,I)*L');
    term2 = -1/(nu+sum(n)+2)*G'*(2*kron(invSig,q) + vec2(q))*G;
    
    fisherinfo.Sigma_ = .5*(term1+term2);
    fisherinfo.n = NaN;
    
    varargout{4} = fisherinfo;
end
%% Optional Parameter Vector Output
if nargout >= 5
    varargout{3} = param;
end
end