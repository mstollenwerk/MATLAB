function [ nLogL, logLcontr, varargout ] = matvsiRiesz2like( Omega_, n, X, varargin )

% Michael Stollenwerk
% michael.stollenwerk@live.com
% 15.03.2022

p = size(Omega_,1);

% [~,iG] = Dmatrix(p);
% I = speye(p);
% L = ELmatrix(p);

Y = matviRiesz2expmat(n);
C_Om = chol(Omega_,'lower');
C_Sig = C_Om/sqrtm(Y);
Sigma_ = C_Sig*C_Sig';
   
[nLogL, logLcontr] = ...
    matviRiesz2like(Sigma_, n, X);
    
if nargout >= 3

    [~, ~, score] = ...
        matviRiesz2like(Sigma_, n, X);

    avg_n = mean(n);
    G = Dmatrix(p);
    invSig = inv(Sigma_);
    fisherinfo_Sigma_iWishart = avg_n/2*G'*kron2(invSig)*G;

    for ii = 1:size(X,3)
        score.rc_paper(:,:,ii) = ...
            ivech(mean(diag(Y))*(fisherinfo_Sigma_iWishart\score.Sigma_(ii,:)'));
    end

    varargout{1} = score;

end

if nargout >= 5

    [~, ~, ~, hessian, param] = ...
        matviRiesz2like(Sigma_, n, X);
    
    varargout{2} = hessian;
    varargout{3} = param;
end

% dOmega_dSigma = iG*kron(C_Sig*Y,I)*L'/(iG*kron(C_Sig,I)*L');
% 
% [nLogL, logLcontr, score, ~, ~, fisherinfo] = ...
%     matviRiesz2like(Sigma_, n, X);
% 
% score.Omega_scaledbyiFish = ...
%     ivech(dOmega_dSigma*(fisherinfo.Sigma_\score.Sigma_'));

end
