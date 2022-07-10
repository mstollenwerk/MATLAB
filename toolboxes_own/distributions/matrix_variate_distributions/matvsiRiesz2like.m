function [ nLogL, logLcontr, varargout ] = ...
    matvsiRiesz2like( Sigma_, nu, X, varargin )
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
    if ~(isempty(Sigma_) && isempty(nu))
        error('Cannot input all_param and any of the parameters individually!')
    end
    all_param = varargin{1};
    Sigma_ = ivechchol(all_param(1 : p_));
    nu = all_param(p_ + 1 : p_ + p);
end
% Checking if Sigma_ is symmetric p.d.
param.Sigma_ = Sigma_;
C = chol(Sigma_,'lower');
param.chol_Sigma_ = vech(C);
param.nu = nu;
diagm = matviRiesz2expmat(nu);
param.all = [param.chol_Sigma_; nu];
%% Log-likelihood computation
logLcontr = NaN(N,1);

term1 = -sum(nu.*log(diag(diagm)))/2;
term2 = -sum(nu)/2*log(2);
term3 = -lgmvgammaln(flipud(nu./2));%-ugmvgammaln(n./2);
term4 = loglpwdet([],nu./2,diag(C)); % upwdet(invS,-n) = lpwdet(S,n)
log_normalizing_constant = term1 + term2 + term3 + term4;

Cr = NaN(p,p,N);
for ii = 1:N
    
    R = X(:,:,ii);
    Cr(:,:,ii) = chol(R,'lower');
    iCz = Cr(:,:,ii)\C;
    
    term5 = loglpwdet([],-(nu+p+1)./2,diag(Cr(:,:,ii))); % upwdet(invS,-n) = lpwdet(S,n)
    term6 = -trace(diagm\(iCz'*iCz))./2;
    
    log_kernel = term5 + term6;

    logLcontr(ii) = log_normalizing_constant + log_kernel;
end
nLogL = -sum(logLcontr);
%% Score computation
if nargout >= 3
    
    score.Sigma_ = NaN(p,p,N);  
    score.nu_originalpdf = NaN(N,p);
    for ii = 1:N

        R = X(:,:,ii);
        
        % General matrix derivative (ignoring symmetry of Sigma_):
        Nabla = - C'\trilHalfDiag(C'*tril(R\C/diagm - C'\diag(nu)))/C; 
        score.SigmaNonSym = Nabla;
        
        % Accounting for symmetry of Sigma_:
        score.Sigma_(:,:,ii) = Nabla+Nabla' - diag(diag(Nabla));

        score.nu_originalpdf(ii,:) = .5*( -log(2) - flip(mvpsi(flip(nu)/2)) ) ...
                - log(diag(Cr(:,:,ii))) + log(diag(C)./sqrt(diag(diagm)));
            
        score.nu_originalpdf_scaled(ii,:) = score.nu_originalpdf(ii,:)' ./ ...
            (.25*flip(mvpsi(flip(nu)/2),1));
                    
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
%     fisherinfo.n = NaN;
%     
%     varargout{4} = fisherinfo;
% end
%% Optional Parameter Vector Output
if nargout >= 5   
    varargout{3} = param;
end
end
