function [ nLogL, logLcontr, varargout ] = ...
    matviRiesz2like( Omega_, nu, X, varargin )
%MATVIRIESZLIKE
%
% Michael Stollenwerk
% michael.stollenwerk@live.com
% 06.09.2021

[p,~,N] = size(X);
p_ = p*(p+1)/2;
narginchk(3,4);
nargoutchk(0,6);
%% Param
if nargin == 4
    if ~(isempty(Omega_) && isempty(nu))
        error('Cannot input all_param and any of the parameters individually!')
    end
    all_param = varargin{1};
    Omega_ = ivechchol(all_param(1 : p_));
    nu = all_param(p_ + 1 : p_ + p);
end
% Checking if Omega_ is symmetric p.d.
param.Omega_ = Omega_;
C = chol(Omega_,'lower');
param.chol_Omega_ = vech(C);
param.n = nu;
param.all = [param.chol_Omega_; nu];
%% Log-likelihood computation
logLcontr = NaN(N,1);

invSig = inv(Omega_);
term1 = -sum(nu)/2*log(2);
term2 = -ugmvgammaln(nu./2);
term3 = loglpwdet([],nu./2,diag(C)); % upwdet(invS,-n) = lpwdet(S,n)
log_normalizing_constant = term1 + term2 + term3;

for ii = 1:N
    
    A = X(:,:,ii);
    
    term4 = loglpwdet(A,-(nu+p+1)./2); % upwdet(invS,-n) = lpwdet(S,n)
    term5 = -trace(Omega_/A)./2;
    
    log_kernel = term4 + term5;

    logLcontr(ii) = log_normalizing_constant + log_kernel;
end
nLogL = -sum(logLcontr);
%% Score computation
if nargout >= 3
    
    score.Omega_ = NaN(N,p_);
    score.n = NaN(N,p);
    
    q = C'\diag(nu)/C;
    
    for ii = 1:N
        invA = inv(X(:,:,ii));
        
        % General matrix derivative (ignoring symmetry of Omega_):
        S = q - invA;
        S = .5*S;
        
        score.Omega_WishFishScaling(:,:,ii) = 2/mean(nu)*Omega_*S*Omega_;

        % Accounting for symmetry of Omega_:
        S = 2*S - diag(diag(S));
        
        % Numerical inaccuracies make S not exactly symmetric
        S = (S+S')./2;
                
        score.Omega_(ii,:) = vech(S);

        % Score nu
        score.n(ii,:) = NaN(p,1);

    end
    
    varargout{1} = score;

end
%% Hessian (Optional Output)

%% Fisher Info (Optional Output)
if nargout >= 6
    
    [G, iG] = Dmatrix(p);
    I = speye(p);
    L = ELmatrix(p);
    
    fishSig = G'*kron(inv(C)',q)*L'/(iG*kron(C,I)*L');
    fishSig = -.5*fishSig;
    fisherinfo.Omega_ = fishSig;
    fisherinfo.n = NaN;
    
    varargout{4} = fisherinfo;
end
%% Optional Parameter Vector Output
if nargout >= 5   
    varargout{3} = param;
end
end
