function [ nLogL, logLcontr, varargout ] = ...
    matvFRieszlike( Omega_, n, nu, X, varargin )
%MATVRIESZLIKE Negative log-likelihood and score of the matrix-variate F distr.
%
% USAGE:
%   [ nLogL, score, param ] = matvFlike( Omega_, n, nu, X )
%   [ nLogL, score, param ] = matvFlike( [], [], [], X, all_param )
%
% INPUTS:
%   X       - Array (p by p). Symmetric p.d. data matrix.
%   Omega_  - Array (p by p). Symmetric p.d. parameter matrix, 
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
    if ~(isempty(Omega_) && isempty(n) && isempty(nu))
        error('Cannot input all_param and any of the parameters individually!')
    end
    all_param = varargin{1};
    Omega_ = ivechchol(all_param(1 : p_));
    n = all_param(p_ + 1 : p_ + p);
    nu = all_param(p_ + p + 1 : p_ + p + p);
end
% Checking if Omega_ is symmetric p.d.
param.Omega_ = Omega_;
COm = chol(Omega_,'lower');
param.chol_Omega_ = vech(COm);
param.n = n;
param.nu = nu;
param.all = [param.chol_Omega_; n; nu];
%% Log-likelihood computation
logLcontr = NaN(N,1);

term1 = ugmvgammaln((n + nu)./2);
term2 = -ugmvgammaln(nu./2);
term3 = -lgmvgammaln(n./2);
term4 = loglpwdet([],nu./2,diag(COm));

log_normalizing_constant = term1 + term2 + term3 + term4;

for ii = 1:N
    term5 = loglpwdet(X(:,:,ii),(n-p-1)./2);
    term6 = loglpwdet(X(:,:,ii) + Omega_, -(n+nu)./2);
    
    log_kernel = term5 + term6;

    logLcontr(ii) = log_normalizing_constant + log_kernel;
end
nLogL = -sum(logLcontr);
%% Score computation
if nargout >= 3
    
    score.Omega_ = NaN(N,p_);
    score.n = NaN(N,p);
    score.nu = NaN(N,p);
    score.n_scaled = NaN(N,p);
    score.nu_scaled = NaN(N,p);
    
    for ii = 1:N
        
        A = X(:,:,ii);
        
        C_Omplusa = chol(Omega_ + A, 'lower');
        S = COm'\diag(nu)/COm - C_Omplusa'\diag(nu+n)/C_Omplusa;
        S = .5*S;
        
        score.Omega_WishFishScaling(:,:,ii) = 2/mean(n)*Omega_*S*Omega_;
        
        % Accounting for symmetry of Omega_:
        S = 2*S - diag(diag(S));
        
        % Numerical inaccuracies make S not exactly symmetric
        S = (S+S')./2;
                
        score.Omega_(ii,:) = vech(S);
        
       
        score.n(ii,:) = .5*flip(mvpsi(flip(n+nu)/2)) - .5*mvpsi(n/2) + log(diag(chol(A,'lower'))) - log(diag(C_Omplusa));
        score.nu(ii,:) = .5*flip(mvpsi(flip(n+nu)/2)) - .5*flip(mvpsi(flip(nu/2))) + log(diag(COm)) - log(diag(C_Omplusa));
        
        score.n_scaled(ii,:) = -score.n(ii,:)'./( .25*flip(mvpsi(flip(n+nu)/2,1)) - .25*mvpsi(n/2,1) );
        score.nu_scaled(ii,:) = -score.nu(ii,:)'./( .25*flip(mvpsi(flip(n+nu)/2,1)) - .25*flip(mvpsi(flip(nu/2),1)) );
        % For scaling see the Fisherinfos below.
        
    end
    
    varargout{1} = score;

end
%% Hessian (Optional Output)
if nargout >= 4
    
    hessian.theta = NaN(2*p,2*p,N);
    for ii = 1:N
        
        hessian.theta(:,:,ii) = ...
            1/4*[diag(flip(mvpsi(flip(n+nu)/2,1)) - mvpsi(n/2,1)), diag(flip(mvpsi(flip(n+nu)/2,1)));
                 diag(flip(mvpsi(flip(n+nu)/2,1))), diag(flip(mvpsi(flip(n+nu)/2,1)) - flip(mvpsi(flip(nu)/2,1)))];
                         
    end
    
    varargout{2} = hessian;
    
end
%% Fisher Info (Optional Output)
if nargout >= 6
   
    %%% I don't have the exact closed form yet, there are some possible
    %%% approximation candidates.
    
    % 1. Approximate with matrix-F fisherinfo, replacing n and nu with the
    % mean(n) and mean(nu). Remember that the FRiesz becomes the matrix-F
    % if all n are equal and all nu are equal. This seems the fairest
    % approximation.
    n_ = mean(n);
    nu_ = mean(nu);    
    c_1 = (n_^2*(nu_-p-2) + 2*n_)/((nu_-p)*(nu_-p-1)*(nu_-p-3));
    c_2 = (n_*(nu_-p-2)+n_^2+n_)/((nu_-p)*(nu_-p-1)*(nu_-p-3));
    c_4 = (n_-p-1)/((n_+nu_-1)*(n_+nu_+2))*((n_-p-2+1/(nu_+n_))*c_2-(1+(n_-p-1)/(n_+nu_))*c_1);
    c_3 = (n_-p-1)/(n_+nu_)*((n_-p-2)*c_2 - c_1)-(n_+nu_+1)*c_4;
    invSigF = inv(Omega_);
    
    G = Dmatrix(p);
    conskron2 = (nu_-(n_+nu_)*(c_3+c_4));
    consvec2 = (n_+nu_)*c_4;
    fishSig = 1/2*G'*(conskron2*kron2(invSigF) - consvec2*vec2(invSigF))*G;
    
     
%     invSig = inv(Omega_);
%     [G, iG] = Dmatrix(p);
%     L = ELmatrix(p);
%     I = speye(p);
%     iC = inv(C);

      % 2. Approximate with part Riesz, part matrix-F fisherinfo. This is
      % motivated by the fact that there are two summands in the FRiesz
      % fisherinfo, the first one we have derived, it is equal to the Riesz
      % fisherinfo, the second one we have not derived its expectation yet.
      % The matrix-F has the same structure, but for the matrix-F we have
      % the expectation of the second summand as well. So we glue the first
      % summand of the FRiezs fisherinfo with the second summand of the
      % matrix-F fisherinfo. This seems a bit aribtrary and does not work
      % well in practice (measured by the fit).
%     n_ = mean(n);
%     nu_ = mean(nu); 
%     c_1 = (n_^2*(nu_-p-2) + 2*n_)/((nu_-p)*(nu_-p-1)*(nu_-p-3));
%     c_2 = (n_*(nu_-p-2)+n_^2+n_)/((nu_-p)*(nu_-p-1)*(nu_-p-3));
%     c_4 = (n_-p-1)/((n_+nu_-1)*(n_+nu_+2))*((n_-p-2+1/(nu_+n_))*c_2-(1+(n_-p-1)/(n_+nu_))*c_1);
%     c_3 = (n_-p-1)/(n_+nu_)*((n_-p-2)*c_2 - c_1)-(n_+nu_+1)*c_4;
% 
%     term_Riesz = 1/2*G'*kron(C'\diag(n),invSig)*L'/(iG*kron(C,I)*L');
%     term_Rest = -(n_ + nu_)/2*G'*((c_3+c_4)*kron2(invSig) + c_4*vec2(invSig))*G;
%     
%     fishSig = term_Riesz + term_Rest;
    
      % 3. Approximate by replacing the expectation of the second summand
      % of the hessian (See 2. Approximation for explanation) hessian,
      % where A is replaced by expectation of A ( so basically do f(E[A]) 
      % instead of E[f(A)] ). This provides good fit, but gives completely
      % out of the range paramter values for the dynamic parameters (~ -3
      % for the score parameter), which are not comparable to other GAS
      % models
%     fishSig = .5*G'*kron(iC',C'\diag(nu)/C)*L'/(G'*kron(C,I)*L')*(G'*G) ...
%               - G'*kron(iC',C'\diag(nu + n)/C)*L'/(G'*kron(C,I)*L')*(G'*G);
          
%     % 4. Approximate by Riesz fisherinfo since for nu->infinity, the
%     FRiesz goes to the Riesz.
%     fishSig = 1/2*G'*kron(C'\diag(n),invSig)*L'/(iG*kron(C,I)*L');
    
%     % 5. Approximate by Wishart fisherinfo since for nu->infinity, all n
%     equal and all nu equal the FRiezs converges to the Wishart. Replace n
%     by mean(n).
%     fishSig = mean(n)/2*G'*kron2(invSig)*G;
    
    fisherinfo.Omega_ = fishSig;
    fisherinfo.dfs = ...
            1/4*[diag(flip(mvpsi(flip(n+nu)/2,1)) - mvpsi(n/2,1)), diag(flip(mvpsi(flip(n+nu)/2,1)));
                 diag(flip(mvpsi(flip(n+nu)/2,1))), diag(flip(mvpsi(flip(n+nu)/2,1)) - flip(mvpsi(flip(nu)/2,1)))];
    
    varargout{4} = fisherinfo;
end
%% Optional Parameter Vector Output
if nargout >= 5
    varargout{3} = param;
end
end
