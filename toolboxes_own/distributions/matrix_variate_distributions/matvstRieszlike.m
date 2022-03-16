function [ nLogL, logLcontr, varargout] = matvstRieszlike( Omega_, n, nu, X, varargin )
% Michael Stollenwerk
% michael.stollenwerk@live.com
% 30.08.2021

Y = diag(n).*nu./(nu - 2);
C_Om = chol(Omega_,'lower');
C_Sig = C_Om/sqrtm(Y);
Sigma_ = C_Sig*C_Sig';
    
[nLogL, logLcontr] = ...
    matvtRieszlike(Sigma_, n, nu, X);
    
if nargout >= 3

    [~, ~, score] = ...
        matvtRieszlike(Sigma_, n, nu, X);

    avg_n = mean(n);
    p = size(Omega_,1);
    G = Dmatrix(p);
    invSig = inv(Sigma_);
    c1 = avg_n/2*(nu+p*avg_n)/(nu+p*avg_n+2);
    c2 = -avg_n^2/2/(nu+p*avg_n+2);
    fisherinfo_Sigma_tWishart = G'*(c1*kron2(invSig) + c2*vec2(invSig))*G;

    for ii = 1:size(X,3)
        score.rc_paper(:,:,ii) = ...
            ivech(mean(diag(Y))*(fisherinfo_Sigma_tWishart\score.Sigma_(ii,:)'));
    end

    varargout{1} = score;

end

if nargout >= 5

    [~, ~, ~, hessian, param] = ...
        matvtRieszlike(Sigma_, n, nu, X);
    
    varargout{2} = hessian;
    varargout{3} = param;
end

%%
% [~,iG] = Dmatrix(p);
% I = speye(p);
% L = ELmatrix(p);
%
% dOmega_dSigma = iG*kron(C_Sig*Y,I)*L'/(iG*kron(C_Sig,I)*L');
% 
% [nLogL, logLcontr, score, ~, ~, fisherinfo] = ...
%     matvtRieszlike(Sigma_, n, nu, X);
% 
% score.Omega_scaledbyiFish = ...
%     ivech(dOmega_dSigma*(fisherinfo.Sigma_\score.Sigma_'));

end
