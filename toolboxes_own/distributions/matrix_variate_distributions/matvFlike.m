function [ nLogL, logLcontr, varargout ] = ...
    matvFlike( Omega_, n, nu, X, varargin )
%MATVFLIKE Negative log-likelihood and score of the matrix-variate F distr.
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
%   df_2    - Double. Second degrees of freedom parameter. 
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
% REFERENCES:
%   [1] Stollenwerk (2022)
%   [2] Gupta and Nagar (2001) - Matrix Variate Distributions, p.167, p.156.
%   [3] Mulder and Pericchi (2018) - The Matrix F Prior for Estimating and
%           Testing Covariance Matrices.
%   [4] Opschoor, Janus, Lucas and Van Dijk (2018) - New HEAVY Models for 
%           Fat-Tailed Realized Covariances and Returns.
%   [5] Bodnar, Okhrin, Parolya and Stollenwerk (2020)
%
% DEPENDENCIES:
%   MVBETALN IVECHCHOL
%
% Michael Stollenwerk
% michael.stollenwerk@live.com
% 16.06.2020

[p,~,N] = size(X);
p_ = p*(p+1)/2;
narginchk(4,5);
nargoutchk(0,6);
%% Param
if nargin == 5
    if ~(isempty(Omega_) && isempty(n) && isempty(nu))
        error('Cannot input all_param and any of the parameters individually!')
    end
    all_param = varargin{1};
    Omega_ = ivechchol(all_param(1 : p_));
    n = all_param(p_ + 1);
    nu = all_param(p_ + 2);
end
% Checking if Omega_ is symmetric p.d.
param.Omega_ = Omega_;
param.chol_Omega_ = vechchol(Omega_);
param.n = n;
param.nu = nu;
param.all = [param.chol_Omega_; n; nu];
%% Log-likelihood computation
logLcontr = NaN(N,1);
term2 = -mvbetaln(n/2, nu/2, p);
term3 = nu/2 * logdet(Omega_);
log_normalizing_constant = term1 + term2 + term3;

for ii = 1:N
    term3 = (n - p - 1)/2*logdet( X(:,:,ii) );
    term4 = -(n + nu)/2*logdet (Omega_ + X(:,:,ii) );
    log_kernel = term3 + term4;

    logLcontr(ii) = log_normalizing_constant + log_kernel;
end
nLogL = -sum(logLcontr);
%% Score computation
if nargout >= 3
    invOmega = inv(Omega_);

    score.Omega_ = NaN(N,p_);
    score.n = NaN(N,1);
    score.nu = NaN(N,1);
    
    for ii = 1:N
        
        S = nu/2*invOmega - (n + nu)/2*inv( Omega_ + X(:,:,ii) );
        
        % Accounting for symmetry of Omega_:
        S = 2*S - diag(diag(S));
        
        % Numerical inaccuracies make S not exactly symmetric
        S = (S+S')./2;
                
        score.Omega_(ii,:) = vech(S);

        % Score nu
        score.df_1(ii) = NaN;
        score.df_2(ii) = NaN;

    end
    
    varargout{1} = score;

end
%% Hessian (Optional Output)
%% Fisher Info (Optional Output)
if nargout >= 6
    
    if nu <= p-3
        error('Fisher Info does not exist for nu <= 3.')
    end
    
    c_1 = (n^2*(nu-p-2) + 2*n)/((nu-p)*(nu-p-1)*(nu-p-3));
    c_2 = (n*(nu-p-2)+n^2+n)/((nu-p)*(nu-p-1)*(nu-p-3));
    c_4 = (n-p-1)/((n+nu-1)*(n+nu+2))*((n-p-2+1/(nu+n))*c_2-(1+(n-p-1)/(n+nu))*c_1);
    c_3 = (n-p-1)/(n+nu)*((n-p-2)*c_2 - c_1)-(n+nu+1)*c_4;
    
    G = Dmatrix(p);
    ckron2 = (nu-(n+nu)*(c_3+c_4));
    cvec2 = (n+nu)*c_4;
    fisherinfo.Omega_ = 1/2*G'*(ckron2*kron2(invOmega) - cvec2*vec2(invOmega))*G;
    fisherinfo.df_1 = NaN;
    fisherinfo.df_2 = NaN;
    
    varargout{4} = fisherinfo;
end
%% Optional Parameter Vector Output
if nargout >= 5  
    varargout{3} = param;
end
end
