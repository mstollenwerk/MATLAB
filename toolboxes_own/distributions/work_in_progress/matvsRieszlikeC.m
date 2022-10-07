function [ nLogL, logLcontr, varargout ] = ...
    matvsRieszlikeC( C, n, Cr, varargin )
%MATVSRIESZLIKE
%
% Michael Stollenwerk
% michael.stollenwerk@live.com
% 20.03.2022

[p,~,N] = size(Cr);
p_ = p*(p+1)/2;
narginchk(3,4);
nargoutchk(0,6);
%% Param
if nargin == 4
    if ~(isempty(C) && isempty(n))
        error('Cannot input all_param and any of the parameters individually!')
    end
    all_param = varargin{1};
    C = ivech(all_param(1 : p_));
    n = all_param(p_ + 1 : p_ + p);
end
% Checking if Sigma_ is symmetric p.d.
param.C = C;
param.n = n;
param.all = [vech(param.C); n];
%% Log-likelihood computation
logLcontr = NaN(N,1);

term1 = sum(n.*log(n))/2;
term2 = -sum(n)/2*log(2);
term3 = -lgmvgammaln(n./2);

for ii = 1:N
    Cz = C(:,:,ii)\Cr(:,:,ii);

    term4 = -loglpwdet([],n./2, diag(C(:,:,ii)));    
    term5 = loglpwdet([],(n-p-1)./2,diag(Cr(:,:,ii)));
    term6 = -trace(diag(n)*(Cz*Cz'))./2;

    logLcontr(ii) = term1 + term2 + term3 + term4 + term5 + term6;
end
nLogL = -sum(logLcontr);
%% Score computation
if nargout >= 3
    
    score.Sigma_ = NaN(p,p,N);  
    score.n_originalpdf = NaN(N,p);
    for ii = 1:N
        
        R = Cr(:,:,ii)*Cr(:,:,ii)';
        C_ = C(:,:,ii);
        
        % General matrix derivative (ignoring symmetry of Sigma_):
        Nabla = C_'\trilHalfDiag(C_'*tril(C_'\diag(n)/C_*R/C_' - C_'\diag(n)))/C_; 
        score.SigmaNonSym = Nabla;
        
        % Accounting for symmetry of Sigma_:
        score.Sigma_(:,:,ii) = Nabla+Nabla' - diag(diag(Nabla));
                
        % Using the chain rule first deriving w.r.t. 
        % Omega_ = C*inv(diag(n))*C' gives the same result, but is much
        % more computationally intensive:
        % ivech(1/2*vec(C'\diag(n)/C*R/C'*diagn/C - C'\diag(n)*diag(n)/C)'*G*dvechCYCt_dvechCCt(C,inv(diag(n))))
        
        score.n_originalpdf(ii,:) = .5*( -log(2) - mvpsi(n/2) ) ...
                        + log(diag(chol(R,'lower'))) - log(diag(C_)./sqrt(n));
        score.n_originalpdf_scaled(ii,:) = score.n_originalpdf(ii,:)' ./ ...
            ( .25*mvpsi(n/2,1) );

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
