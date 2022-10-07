function [ nLogL, logLcontr, varargout ] = ...
    matvsiFRiesz2like( Sigma_, n, nu, R, varargin )

% Michael Stollenwerk
% michael.stollenwerk@live.com
% 28.04.2022

[p,~,N] = size(R);
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
    n = all_param(p_ + 1 : p_ + p);
    nu = all_param(p_ + p + 1 : p_ + p + p);
end
% Checking if Sigma_ is symmetric p.d.
param.Sigma_ = Sigma_;
C = chol(Sigma_,'lower');
A = matviFRiesz2expmat(n,nu);
C_Omega = C/sqrt(A);
iC_Omega = inv(C_Omega);
param.chol_Sigma_ = vechchol(Sigma_);
param.nu = nu;
param.n = n;
param.all = [param.chol_Sigma_; n; nu];
%% Log-likelihood computation
logLcontr = NaN(N,1);

term1 = lgmvgammaln((nu + n)./2);
term2 = -lgmvgammaln(n./2);
term3 = -lgmvgammaln(flip(nu)./2);
term4 = loglpwdet([],-n./2,diag(C_Omega));

log_normalizing_constant = term1 + term2 + term3 + term4;

C_B = NaN(p,p,N);
C_R = NaN(p,p,N);
for ii = 1:N
    
    C_R(:,:,ii) = chol(R(:,:,ii),'lower');
    iCR = inv(C_R(:,:,ii));
    
    B = inv(iCR'*iCR + iC_Omega'*iC_Omega);
    C_B(:,:,ii) = chol(B,'lower');
    
    term5 = loglpwdet([],-(nu+p+1)./2, diag(C_R(:,:,ii)));
    term6 = loglpwdet([], (n+nu)./2, diag(C_B(:,:,ii)));
    
    log_kernel = term5 + term6;

    logLcontr(ii) = log_normalizing_constant + log_kernel;
end
nLogL = -sum(logLcontr);
%% Optional Parameter Vector Output
if nargout >= 3
    
    score.Sigma_ = NaN(p,p,N);
    score.n = NaN(N,p);
    score.nu = NaN(N,p);
    
    for ii = 1:N
        
        Nabla = -(C'\trilHalfDiag(C'*tril(C'\diag(n) - C'\A/C*C_B(:,:,ii)*diag(n+nu)*C_B(:,:,ii)'/C'))/C);
        score.SigmaNonSym(:,:,ii) = Nabla;

        % Accounting for symmetry of Sigma_:    
        score.Sigma_(:,:,ii) = Nabla+Nabla' - diag(diag(Nabla));
                
        score.n_originalpdf(ii,:) = .5*(mvpsi((n+nu)/2) - mvpsi(n/2)) ...
                      - log(diag(C_Omega)) + log(diag(C_B(:,:,ii)));
        score.n_originalpdf_scaled(ii,:) = -score.n_originalpdf(ii,:)'./( .25*( mvpsi((n+nu)/2,1) - mvpsi(n/2,1) ) ); 
        score.nu_originalpdf(ii,:) = .5*(mvpsi((n+nu)/2) - flip(mvpsi(flip(nu)/2))) ...
                      - log(diag(C_R(:,:,ii))) + log(diag(C_B(:,:,ii)));
        score.nu_originalpdf_scaled(ii,:) = -score.nu_originalpdf(ii,:)'./( .25*(mvpsi((n+nu)/2,1) - flip(mvpsi(flip(nu)/2,1))) );

    end
    
    varargout{1} = score;
    
end
if nargout >= 4
    varargout{2} = NaN;
end
if nargout >= 5
    varargout{3} = param;
end
end
