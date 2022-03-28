function [ nLogL, logLcontr, varargout ] = ...
    matvsRieszlike( Sigma_, n, X, varargin )
%MATVSRIESZLIKE
%
% Michael Stollenwerk
% michael.stollenwerk@live.com
% 20.03.2022

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

term1 = sum(n.*log(n))/2;
term2 = -sum(n)/2*log(2);
term3 = -lgmvgammaln(n./2);
term4 = -loglpwdet([],n./2, diag(C));

log_normalizing_constant = term1 + term2 + term3 + term4;

for ii = 1:N
    R = X(:,:,ii);
    Cr = chol(R,'lower');
    Cz = C\Cr;
   
    term5 = loglpwdet([],(n-p-1)./2,diag(Cr));
    term6 = -trace(diag(n)*(Cz*Cz'))./2;
    
    log_kernel = term5 + term6;

    logLcontr(ii) = log_normalizing_constant + log_kernel;
end
nLogL = -sum(logLcontr);
%% Score computation
if nargout >= 3
    
    score.Sigma_ = NaN(p,p,N);    
    for ii = 1:N
        
        R = X(:,:,ii);
        
        % General matrix derivative (ignoring symmetry of Sigma_):
        Nabla = C'\trilHalfDiag(C'*tril(C'\diag(n)/C*R/C' - C'\diag(n)))/C; 
        score.SigmaNonSym = Nabla;
        
        % Accounting for symmetry of Sigma_:
        score.Sigma_(:,:,ii) = Nabla+Nabla' - diag(diag(Nabla));
                
        % Using the chain rule first deriving w.r.t. 
        % Omega_ = C*inv(diag(n))*C' gives the same result, but is much
        % more computationally intensive:
        % ivech(1/2*vec(C'\diag(n)/C*R/C'*diagn/C - C'\diag(n)*diag(n)/C)'*G*dvechCYCt_dvechCCt(C,inv(diag(n))))

    end
    
    varargout{1} = score;

end
%% Hessian (Optional Output)
%% Fisher Info (Optional Output)
% if nargout >= 6
%     
%     [G, iG] = Dmatrix(p);
%     I = speye(p);
%     L = ELmatrix(p);
%     
%     invSig = inv(Sigma_);
%     
%     fisherinfo.Sigma_ = NaN;
%     fisherinfo.n = NaN;
%     
%     varargout{4} = fisherinfo;
% end
%% Optional Parameter Vector Output
if nargout >= 5
    varargout{3} = param;
end
end
