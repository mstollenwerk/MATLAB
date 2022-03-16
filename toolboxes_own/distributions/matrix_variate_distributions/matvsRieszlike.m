function [ nLogL, logLcontr, varargout ] = matvsRieszlike( Omega_, n, X, varargin )
%MATVWISHLIKE
%
% USAGE:
%   
%
% INPUTS:
%   
%
% OUTPUTS:
%   
%  See also 
%
% COMMENTS:
%   
% REFERENCES:
%      [1]                    
%
% DEPENDENCIES:
%
%
% Michael Stollenwerk
% michael.stollenwerk@live.com
% 30.08.2021

Y = diag(n);
C_Om = chol(Omega_,'lower');
C_Sig = C_Om/sqrtm(Y);
Sigma_ = C_Sig*C_Sig';

[nLogL, logLcontr] = ...
    matvRieszlike(Sigma_, n, X);
    
if nargout >= 3
    
    [~, ~, score] = ...
        matvRieszlike(Sigma_, n, X);

    avg_n = mean(n);
    p = size(Omega_,1);
    G = Dmatrix(p);
    invSig = inv(Sigma_);
    fisherinfo_Sigma_Wishart = avg_n/2*G'*kron(invSig,invSig)*G;
    
    for ii = 1:size(X,3)
        score.rc_paper(:,:,ii) = ...
            ivech(avg_n*(fisherinfo_Sigma_Wishart\score.Sigma_(ii,:)'));
    end

    varargout{1} = score;

end

if nargout >= 5

    [~, ~, ~, hessian, param] = ...
        matvRieszlike(Sigma_, n, X);
    
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
%     matvRieszlike(Sigma_, n, X);
% 
% score.Omega_scaledbyiFish = ...
%     ivech(dOmega_dSigma*(fisherinfo.Sigma_\score.Sigma_'));

end
