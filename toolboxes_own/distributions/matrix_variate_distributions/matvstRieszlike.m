function [ nLogL, logLcontr, varargout ] = ...
    matvstRieszlike( Sigma_, n, nu, X, varargin )
%MATVSRIESZTLIKE
%
% Michael Stollenwerk
% michael.stollenwerk@live.com
% 11.02.2021

[p,~,N] = size(X);
p_ = p*(p+1)/2;
narginchk(4,5); %%%%%%%
nargoutchk(0,6);
%% Param
if nargin >= 5 %%%%%%%
    if ~(isempty(Sigma_) && isempty(nu) && isempty(n))
        error('Cannot input all_param and any of the parameters individually!')
    end
    all_param = varargin{1};
    Sigma_ = ivechchol(all_param(1 : p_));
    n = all_param(p_ + 1 : p_ + p);
    nu = all_param(p_ + p + 1);
end
% Checking if Sigma_ is symmetric p.d.
param.Sigma_ = Sigma_;
C = chol(Sigma_,'lower');
param.chol_Sigma_ = vech(C);
param.nu = nu;
param.n = n;
param.all = [param.chol_Sigma_; n; nu];
%% Log-likelihood computation
logLcontr = NaN(N,1);

term1 = sum(n.*log(n))/2;
term2 = gammaln((nu + sum(n))/2);
term3 = -gammaln(nu/2);
term4 = -lgmvgammaln(n./2);
term5 = -sum(n)/2*log(nu-2);
term6 = -loglpwdet([],n./2,diag(C));

log_normalizing_constant = term1 + term2 + term3 + term4 + term5 + term6;

for ii = 1:N
    R = X(:,:,ii);
    Cr = chol(R,'lower');
    Cz = C\Cr;
    
    term7 = loglpwdet([],(n-p-1)./2,diag(Cr));
    term8 = -(nu + sum(n))/2*log(1 + trace(diag(n)*(Cz*Cz'))/(nu-2));
    
    log_kernel = term7 + term8;

    logLcontr(ii) = log_normalizing_constant + log_kernel;
end
nLogL = -sum(logLcontr);
%% Score computation
if nargout >= 3
        
    score.Sigma_ = NaN(p,p,N);    
    for ii = 1:N
        
        R = X(:,:,ii);
        diagnZ = diag(n)/C*R/C';
        % General matrix derivative (ignoring symmetry of Sigma_):
        Nabla = C'\trilHalfDiag(C'*tril((nu+sum(n))/(nu-2+trace(diagnZ))*(C'\diagnZ) - C'\diag(n)))/C; 
        
        % Accounting for symmetry of Sigma_:
        score.Sigma_(:,:,ii) = Nabla+Nabla' - diag(diag(Nabla));

    end
    
    varargout{1} = score;

end
%% Hessian (Optional Output)

%% Fisher Info (Optional Output)
if nargout >= 6
%     
%     G = Dmatrix(p);
%     I = speye(p);
%     L = ELmatrix(p);
%     fisherinfo.n = NaN;
%     
%     varargout{4} = fisherinfo;
end
%% Optional Parameter Vector Output
if nargout >= 5   
    varargout{3} = param;
end
end