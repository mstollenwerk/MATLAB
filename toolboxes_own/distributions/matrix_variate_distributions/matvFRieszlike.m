function [ nLogL, logLcontr, varargout ] = ...
    matvFRieszlike( Sigma_, n, nu, X, varargin )
%MATVRIESZLIKE Negative log-likelihood and score of the matrix-variate F distr.
%
% USAGE:
%   [ nLogL, score, param ] = matvFlike( Sigma_, n, nu, X )
%   [ nLogL, score, param ] = matvFlike( [], [], [], X, all_param )
%
% INPUTS:
%   X       - Array (p by p). Symmetric p.d. data matrix.
%   SIGMA_  - Array (p by p). Symmetric p.d. parameter matrix, 
%             regulates covariance. 
%   n    - Double. First degrees of freedom parameter. 
%   nu    - Double. Second degrees of freedom parameter. 
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
    if ~(isempty(Sigma_) && isempty(n) && isempty(nu))
        error('Cannot input all_param and any of the parameters individually!')
    end
    all_param = varargin{1};
    Sigma_ = ivechchol(all_param(1 : p_));
    n = all_param(p_ + 1 : p_ + p);
    nu = all_param(p_ + p + 1 : p_ + p + p);
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

term1 = ugmvgammaln((n + nu)./2);
term2 = -ugmvgammaln(nu./2);
term3 = -lgmvgammaln(n./2);
term4 = loglpwdet([],nu./2,diag(C));

log_normalizing_constant = term1 + term2 + term3 + term4;

for ii = 1:N
    term5 = loglpwdet(X(:,:,ii),(n-p-1)./2);
    term6 = loglpwdet(X(:,:,ii) + Sigma_, -(n+nu)./2);
    
    log_kernel = term5 + term6;

    logLcontr(ii) = log_normalizing_constant + log_kernel;
end
nLogL = -sum(logLcontr);
%% Score computation
if nargout >= 3
    
    score.Sigma_ = NaN(N,p_);
    score.n = NaN(N,p);
    score.nu = NaN(N,p);
    
    for ii = 1:N
        
        A = X(:,:,ii);
        
        C_sigplusa = chol(Sigma_ + A, 'lower');
        S = C'\diag(nu)/C - C_sigplusa'\diag(nu+n)/C_sigplusa;
        S = .5*S;
        
        score.Sigma_WishFishScaling(:,:,ii) = 2/mean(n)*Sigma_*S*Sigma_;
        
        % Accounting for symmetry of Sigma_:
        S = 2*S - diag(diag(S));
        
        % Numerical inaccuracies make S not exactly symmetric
        S = (S+S')./2;
                
        score.Sigma_(ii,:) = vech(S);

    end
    
    varargout{1} = score;

end
%% Hessian (Optional Output)
%% Fisher Info (Optional Output)
if nargout >= 6
    
    %%% THIS IS ONLY AN APPROXIMATION TO THE FISHER INFO WRT SIGMA!!! %%%
%     n_ = mean(n);
%     nu_ = mean(nu);    
%     c_1 = (n_^2*(nu_-p-2) + 2*n_)/((nu_-p)*(nu_-p-1)*(nu_-p-3));
%     c_2 = (n_*(nu_-p-2)+n_^2+n_)/((nu_-p)*(nu_-p-1)*(nu_-p-3));
%     c_4 = (n_-p-1)/((n_+nu_-1)*(n_+nu_+2))*((n_-p-2+1/(nu_+n_))*c_2-(1+(n_-p-1)/(n_+nu_))*c_1);
%     c_3 = (n_-p-1)/(n_+nu_)*((n_-p-2)*c_2 - c_1)-(n_+nu_+1)*c_4;
%     invSigF = inv(Sigma_);
%     
%     G = Dmatrix(p);
%     ckron2 = (nu_-(n_+nu_)*(c_3+c_4));
%     cvec2 = (n_+nu_)*c_4;
%     fisherinfo.Sigma_ = 1/2*G'*(ckron2*kron2(invSigF) - cvec2*vec2(invSigF))*G;
    G = Dmatrix(p);
    L = ELmatrix(p);
    I = speye(p);
    
    iC = inv(C);
    
    %%% THIS IS ONLY AN APPROXIMATION TO THE FISHER INFO WRT SIGMA!!! %%%
    fishSig = .5*G'*kron(iC',C'\diag(nu)/C)*L'/(G'*kron(C,I)*L')*(G'*G) ...
              - G'*kron(iC',C'\diag(nu + n)/C)*L'/(G'*kron(C,I)*L')*(G'*G);
    
    fisherinfo.Sigma_ = fishSig;
    fisherinfo.df_1 = NaN;
    fisherinfo.df_2 = NaN;
    
    varargout{4} = fisherinfo;
end
%% Optional Parameter Vector Output
if nargout >= 5
    varargout{3} = param;
end
end
