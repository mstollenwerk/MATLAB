function [ nLogL, logLcontr, varargout ] = ...
    matvsFlike( Sigma_, n, nu, X, varargin )
%MATVSFLIKE Negative log-likelihood of the standardized matrix-variate F
%distribution.
%
% USAGE:
%   [ nLogL, score, param ] = matvFlike( Sigma_, df_1, df_2, X )
%   [ nLogL, score, param ] = matvFlike( [], [], [], X, all_param )
%
% INPUTS:
%   X       - Array (p by p). Symmetric p.d. data matrix.
%   SIGMA_  - Array (p by p). Symmetric p.d. parameter matrix, 
%             expected value matrix. 
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
    if ~(isempty(Sigma_) && isempty(n) && isempty(nu))
        error('Cannot input all_param and any of the parameters individually!')
    end
    all_param = varargin{1};
    Sigma_ = ivechchol(all_param(1 : p_));
    n = all_param(p_ + 1);
    nu = all_param(p_ + 2);
end
% Checking if Sigma_ is symmetric p.d.
param.Sigma_ = Sigma_;
C = chol(Sigma_,'lower');
param.chol_Sigma_ = vech(C);
param.n = n;
param.nu = nu;
param.all = [param.chol_Sigma_; n; nu];
%% Log-likelihood computation
logLcontr = NaN(N,1);
term1 =   p*n/2*log(n/(nu-p-1));
term2 = - mvbetaln(n/2, nu/2, p);
term3 = - n/2*logdet(Sigma_);
log_normalizing_constant = term1 + term2 + term3;

I = eye(p);
for ii = 1:N
    R = X(:,:,ii);
    
    term4 = (n - p - 1)/2*logdet(R);
    term5 = -(n + nu)/2*log(det( I + n/(nu-p-1)*(Sigma_\R) ));
    log_kernel = term4 + term5;

    logLcontr(ii) = log_normalizing_constant + log_kernel;
end
nLogL = -sum(logLcontr);
%% Score computation
if nargout >= 3
    invSig = inv(Sigma_);

    score.Sigma_ = NaN(p,p,N);
    score.n_originalpdf = NaN(N,1);
    score.nu_originalpdf = NaN(N,1);
    
    for ii = 1:N
        
        B = Sigma_ + X(:,:,ii)*n/(nu-p-1);
        S = nu*invSig - (n + nu)*inv( B );
        S = .5.*S;
        
        score.SigmaNonSym = S;
        
        % Accounting for symmetry of Sigma_:
        S = S+S' - diag(diag(S));
           
        score.Sigma_(:,:,ii) = S;
        
        Omega_ = matvStandardize('F',Sigma_,[n,nu]);
        score.n_originalpdf(ii) = .5*( sum(mvpsi(ones(p,1)*(n+nu)/2)) - sum(mvpsi(ones(p,1)*n/2)) ...
                           + logdet(X(:,:,ii)) - logdet(Omega_ + X(:,:,ii)) );
        score.nu_originalpdf(ii) = .5*( sum(mvpsi(ones(p,1)*(n+nu)/2)) - sum(mvpsi(ones(p,1)*nu/2)) ...
                            + logdet(Omega_) - logdet(Omega_ + X(:,:,ii)) );
                        
        fisherinfo.n_originalpdf(ii) = -.25*( sum(mvpsi(ones(p,1)*(n+nu)/2,1)) - sum(mvpsi(ones(p,1)*n/2,1)) );
        fisherinfo.nu_originalpdf(ii) = -.25*( sum(mvpsi(ones(p,1)*(n+nu)/2,1)) - sum(mvpsi(ones(p,1)*nu/2,1)) );
        
        score.n_originalpdf_scaled(ii) = score.n_originalpdf(ii)/fisherinfo.n_originalpdf(ii);
        score.nu_originalpdf_scaled(ii) = score.nu_originalpdf(ii)/fisherinfo.nu_originalpdf(ii);                        

    end
    
    varargout{1} = score;

end
%% Hessian (Optional Output)
%% Fisher Info (Optional Output)
if nargout >= 6
    
    c_1 = (n^2*(nu-p-2) + 2*n)/((nu-p)*(nu-p-1)*(nu-p-3));
    c_2 = (n*(nu-p-2)+n^2+n)/((nu-p)*(nu-p-1)*(nu-p-3));
    c_4 = (n-p-1)/((n+nu-1)*(n+nu+2))*((n-p-2+1/(nu+n))*c_2-(1+(n-p-1)/(n+nu))*c_1);
    c_3 = (n-p-1)/(n+nu)*((n-p-2)*c_2 - c_1)-(n+nu+1)*c_4;
    
    G = Dmatrix(p);
    ckron2 = (nu-(n+nu)*(c_3+c_4));
    cvec2 = (n+nu)*c_4;
    fisherinfo.Sigma_ = 1/2*G'*(ckron2*kron2(invSig) - cvec2*vec2(invSig))*G;
    fisherinfo.df_1 = NaN;
    fisherinfo.df_2 = NaN;
    
    varargout{4} = fisherinfo;
end
%% Optional Parameter Vector Output
if nargout >= 5  
    varargout{3} = param;
end
end
