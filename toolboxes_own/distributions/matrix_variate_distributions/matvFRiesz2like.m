function [ nLogL, logLcontr, varargout ] = ...
    matvFRiesz2like( Sigma_, n, nu, X, varargin )

% Michael Stollenwerk
% michael.stollenwerk@live.com
% 31.08.2021

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
U = cholU(Sigma_);
param.chol_Sigma_ = vechchol(Sigma_);
param.n = n;
param.nu = nu;
param.all = [param.chol_Sigma_; n; nu];
%% Log-likelihood computation
logLcontr = NaN(N,1);

term1 = lgmvgammaln((n + nu)./2);
term2 = -lgmvgammaln(nu./2);
term3 = -ugmvgammaln(n./2);
term4 = logupwdet([],nu./2,diag(U));

log_normalizing_constant = term1 + term2 + term3 + term4;

for ii = 1:N
    term5 = logupwdet(X(:,:,ii),(n-p-1)./2);
    term6 = logupwdet(X(:,:,ii) + Sigma_, -(n+nu)./2);
    
    log_kernel = term5 + term6;

    logLcontr(ii) = log_normalizing_constant + log_kernel;
end
nLogL = -sum(logLcontr);
%% Score computation
if nargout >= 3
    error('There might be a mistake either in score or in fisher info')
    score.Sigma_ = NaN(N,p_);
    score.n = NaN(N,p);
    score.nu = NaN(N,p);
    
    for ii = 1:N
        
        A = X(:,:,ii);
        
        U_sigplusa = cholU(Sigma_ + A);
        S = U'\diag(nu)/U - U_sigplusa'\diag(nu+n)/U_sigplusa;
        S = .5*S;
        
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
    error('There might be a mistake either in score or in fisher info')
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
    
    iC = inv(U);
    
    %%% THIS IS ONLY AN APPROXIMATION TO THE FISHER INFO WRT SIGMA!!! %%%
    fishSig = .5*G'*kron(iC',U'\diag(nu)/U)*L'/(G'*kron(U,I)*L')*(G'*G) ...
              - G'*kron(iC',U'\diag(nu + n)/U)*L'/(G'*kron(U,I)*L')*(G'*G);
    
    fisherinfo.Sigma_ = fishSig;
    fisherinfo.df_1 = NaN;
    fisherinfo.df_2 = NaN;
    
    varargout{4} = fisherinfo;
end
%% Optional Parameter Vector Output
if nargout == 5
    varargout{3} = param;
end
end
