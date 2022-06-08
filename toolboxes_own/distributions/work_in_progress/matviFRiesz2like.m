function [ nLogL, logLcontr, varargout ] = ...
    matviFRiesz2like( Omega_, n, nu, X, varargin )

% Michael Stollenwerk
% michael.stollenwerk@live.com
% 28.04.2022

[p,~,N] = size(X);
p_ = p*(p+1)/2;
narginchk(4,5); %%%%%%%
nargoutchk(0,6);
%% Param
if nargin == 5 %%%%%%%
    if ~(isempty(Omega_) && isempty(nu) && isempty(n))
        error('Cannot input all_param and any of the parameters individually!')
    end
    all_param = varargin{1};
    Omega_ = ivechchol(all_param(1 : p_));
    n = all_param(p_ + 1 : p_ + p);
    nu = all_param(p_ + p + 1 : p_ + p + p);
end
% Checking if Sigma_ is symmetric p.d.
param.Omega_ = Omega_;
C_Omega = chol(Omega_,'lower');
param.chol_Omega_ = vechchol(Omega_);
param.nu = nu;
param.n = n;
param.all = [param.chol_Omega_; n; nu];
%% Log-likelihood computation
logLcontr = NaN(N,1);

term1 = lgmvgammaln((nu + n)./2);
term2 = -lgmvgammaln(n./2);
term3 = -lgmvgammaln(flip(nu)./2);
term4 = loglpwdet([],-n./2,diag(C_Omega));

log_normalizing_constant = term1 + term2 + term3 + term4;

for ii = 1:N
    term5 = loglpwdet(X(:,:,ii),-(nu+p+1)./2);
    term6 = logupwdet(inv(X(:,:,ii)) + inv(Omega_), -(n+nu)./2);
    
    log_kernel = term5 + term6;

    logLcontr(ii) = log_normalizing_constant + log_kernel;
end
nLogL = -sum(logLcontr);
%% Optional Parameter Vector Output
if nargout >= 3
    
    score.Sigma_ = NaN(p,p,N);
    
    for ii = 1:N
        
        B = inv(inv(X(:,:,ii)) + inv(Omega_));
        C_B = chol(B,'lower');
        
        % General matrix derivative (ignoring symmetry of Sigma_):
        S = -1/2*(C_Omega'\diag(n)/C_Omega - Omega_\C_B*diag(n+nu)*C_B'/Omega_);
        
%         score.Sigma_WishFishScaling(:,:,ii) = 2/mean(n)*Sigma_*S*Sigma_;
        
        % Accounting for symmetry of Sigma_:
        S = 2*S - diag(diag(S));
        
        % Numerical inaccuracies make S not exactly symmetric
        S = (S+S')./2;
                
        score.Omega_(:,:,ii) = S;

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
