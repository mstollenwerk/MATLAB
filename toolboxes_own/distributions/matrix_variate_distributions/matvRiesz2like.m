function [ nLogL, logLcontr, varargout ] = ...
    matvRiesz2like( Sigma_, n, X, varargin )

% Michael Stollenwerk
% michael.stollenwerk@live.com
% 31.08.2021

[p,~,N] = size(X);
p_ = p*(p+1)/2;
narginchk(3,4);
nargoutchk(0,6);
%% Param
if nargin == 4
    if ~(isempty(Sigma_) && isempty(n))
        error('Cannot input all_param and any of the parameters individually!')
    end
    all_param = varargin{1};
    Sigma_ = ivechchol(all_param(1 : p_));
    n = all_param(p_ + 1 : p_ + p);
end
% Checking if Sigma_ is symmetric p.d.
param.Sigma_ = Sigma_;
param.chol_Sigma_ = vechchol(Sigma_);
param.n = n;
param.all = [param.chol_Sigma_; n];
%% Log-likelihood computation
logLcontr = NaN(N,1);

Cdot = chol(inv(Sigma_),'lower');
iCdot = inv(Cdot);

term1 = -sum(n)/2*log(2);
term2 = -ugmvgammaln(n./2);
term3 = logupwdet([],-n./2, diag(iCdot));

log_normalizing_constant = term1 + term2 + term3;

for ii = 1:N
    
    A = X(:,:,ii);
    
    term4 = logupwdet(A,(n-p-1)./2);
    term5 = -trace(Sigma_\A)./2;
    
    log_kernel = term4 + term5;

    logLcontr(ii) = log_normalizing_constant + log_kernel;
end
nLogL = -sum(logLcontr);
%% Score computation
if nargout >= 3
    
    score.Sigma_ = NaN(N,p_);
    score.n = NaN(N,p);
    
    for ii = 1:N
        
        A = X(:,:,ii);
        
        % General matrix derivative (ignoring symmetry of Sigma_):
        S = 1/2*(Sigma_\A/Sigma_ - Cdot*diag(n)*Cdot');
        
        score.Sigma_WishFishScaling(:,:,ii) = 2/mean(n)*Sigma_*S*Sigma_;
        
        % Accounting for symmetry of Sigma_:
        S = 2*S - diag(diag(S));
        
        % Numerical inaccuracies make S not exactly symmetric
        S = (S+S')./2;
                
        score.Sigma_(ii,:) = vech(S);

        % Score nu
        score.n(ii,:) = NaN(p,1);

    end
    
    varargout{1} = score;

end
%% Hessian (Optional Output)
%% Fisher Info (Optional Output)
if nargout >= 6
    
    [G, iG] = Dmatrix(p);
    L = ELmatrix(p);

    
    fisherinfo.Sigma_ = ...
        1/2*G'*kron(Cdot,Cdot*diag(n)*iCdot)*L'/(iG*kron(iCdot',Sigma_)*L');
    fisherinfo.n = NaN;
    
    varargout{4} = fisherinfo;
end
%% Optional Parameter Vector Output
if nargout >= 5
    varargout{3} = param;
end
end
