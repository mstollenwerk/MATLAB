function [ nLogL, logLcontr, varargout ] = ...
    matvsitRiesz2like( Sigma_, n, nu, X, varargin )
%MATVITRIESZ2LIKE
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
    if ~(isempty(Sigma_) && isempty(n) && isempty(nu))
        error('Cannot input all_param and any of the parameters individually!')
    end
    all_param = varargin{1};
    Sigma_ = ivechchol(all_param(1 : p_));
    n = all_param(p_ + 1);    
    nu = all_param(p_ +  2 : p_ + 1 + p);
end
% Checking if Sigma_ is symmetric p.d.
param.Sigma_ = Sigma_;
C = chol(Sigma_,'lower');
param.chol_Sigma_ = vech(C);
param.n = n;
param.nu = nu;
diagm = matviRiesz2expmat(nu);
param.all = [param.chol_Sigma_; n; nu];
%% Log-likelihood computation
logLcontr = NaN(N,1);

term1 = -sum(nu.*log(diag(diagm)))/2;
term2 = gammaln((n + sum(nu))/2);
term3 = -gammaln(n/2);
term4 = -ugmvgammaln(nu./2);
term5 = -sum(nu)/2*log(n);
term6 = loglpwdet([],nu./2,diag(C)); % upwdet(invS,-n) = lpwdet(S,n)

log_normalizing_constant = term1 + term2 + term3 + term4 + term5 + term6;

for ii = 1:N
    
    R = X(:,:,ii);
    Cr = chol(R,'lower');
    iCz = Cr\C;
    
    term7 = loglpwdet([],-(nu+p+1)./2,diag(Cr)); % upwdet(invS,-n) = lpwdet(S,n)
    term8 = -(n + sum(nu))/2*log(1 + trace(diagm\(iCz'*iCz))/n);
    
    log_kernel = term7 + term8;

    logLcontr(ii) = log_normalizing_constant + log_kernel;
end
nLogL = -sum(logLcontr);
%% Score computation
if nargout >= 3
    
    score.Sigma_ = NaN(p,p,N);
    score.n_originalpdf = NaN(N,1);
    score.nu_originalpdf = NaN(N,p);
    for ii = 1:N

        R = X(:,:,ii);
        
        Om = C/diagm*C';
        COm = chol(Om,'lower');
        trOmInvR = trace(Om/R);
        % General matrix derivative (ignoring symmetry of Sigma_):
        Nabla = - C'\trilHalfDiag(C'*tril( ...
            (n + sum(nu))/(n+trOmInvR)*inv(R)*C*inv(diagm) - C'\diag(nu) ...
            ))/C; 
        score.SigmaNonSym = Nabla;
        
        % Accounting for symmetry of Sigma_:
        score.Sigma_(:,:,ii) = Nabla+Nabla' - diag(diag(Nabla));

        score.n_originalpdf(ii) = .5*( psi((n+sum(nu))/2) - psi(n/2) - sum(nu)/n - ...
                           log(1+trOmInvR/n) + (n+sum(nu))*(trOmInvR/n^2)/(1+trOmInvR/n) );
        score.nu_originalpdf(ii,:) = .5*( psi((n+sum(nu))/2) - flip(mvpsi(flip(nu)/2)) - log(n+trOmInvR) ) ...
                         + log(diag(COm)) - log(diag(chol(R,'lower'))); 
        score.n_originalpdf_scaled(ii) = -score.n_originalpdf(ii) ./ ...
            ( .25*psi(1,(n+sum(nu))/2) - .25*psi(1,n/2) + .5*(sum(nu)+n+4)/(sum(nu)+n+2)*sum(nu)/(sum(nu)+n)/n );
        score.nu_originalpdf_scaled(ii,:) = -score.nu_originalpdf(ii,:)' ./ ...
            ( .25*( psi(1,(n+sum(nu))/2) - flip(mvpsi(flip(nu)/2,1)) ) );
        
    end
    
    varargout{1} = score;

end
%% Hessian (Optional Output)

%% Fisher Info (Optional Output)
% if nargout >= 6
%     
%     [G,iG] = Dmatrix(p);
%     I = speye(p);
%     L = ELmatrix(p);
%     
%     varargout{4} = fisherinfo;
% end
%% Optional Parameter Vector Output
if nargout >= 5
    varargout{3} = param;
end
end