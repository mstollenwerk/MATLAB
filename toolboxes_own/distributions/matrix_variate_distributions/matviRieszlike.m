function [ nLogL, logLcontr, varargout ] = ...
    matviRieszlike( Sigma_, n, X, varargin )
%MATVIRIESZLIKE
%
% Michael Stollenwerk
% michael.stollenwerk@live.com
% 11.02.2021

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
C = chol(Sigma_,'lower');
param.chol_Sigma_ = vech(C);
param.n = n;
param.all = [param.chol_Sigma_; n];
%% Log-likelihood computation
logLcontr = NaN(N,1);

invSig = inv(Sigma_);
Cdot = chol(invSig,'lower');
term1 = -sum(n)/2*log(2);
term2 = -lgmvgammaln(n./2);
term3 = -loglpwdet([],n./2,diag(Cdot));
log_normalizing_constant = term1 + term2 + term3;

for ii = 1:N
    
    term4 = loglpwdet(inv(X(:,:,ii)),(n+p+1)./2);
    term5 = -trace(Sigma_/X(:,:,ii))./2;
    
    log_kernel = term4 + term5;

    logLcontr(ii) = log_normalizing_constant + log_kernel;
end
nLogL = -sum(logLcontr);
%% Score computation
if nargout >= 3
    
    score.Sigma_ = NaN(N,p_);
    score.n = NaN(N,p);
    
    for ii = 1:N
        invA = inv(X(:,:,ii));
        
        % General matrix derivative (ignoring symmetry of Sigma_):
        S = Cdot*diag(n)*Cdot' - invA;
        S = .5*S;

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
    
    G = Dmatrix(p);
    I = speye(p);
    L = ELmatrix(p);
    
    fishSig = G'*kron(Cdot*diag(n),I)*L'/(G'*kron(Cdot,I)*L')*G'*kron2(invSig)*G;
    fishSig = -.5*fishSig;
    fisherinfo.Sigma_ = fishSig;
    fisherinfo.n = NaN;
    
    varargout{4} = fisherinfo;
end
%% Optional Parameter Vector Output
if nargout >= 5   
    varargout{3} = param;
end
end
