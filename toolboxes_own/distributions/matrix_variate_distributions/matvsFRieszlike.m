function [ nLogL, logLcontr, varargout ] = matvsFRieszlike( Omega_, n, nu, X, varargin )
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

Y = matvFRieszexpmat(n,nu);
C_Om = chol(Omega_,'lower');
C_Sig = C_Om/sqrtm(Y);
Sigma_ = C_Sig*C_Sig';
    
[nLogL, logLcontr] = ...
    matvFRieszlike(Sigma_, n, nu, X);
    
if nargout >= 3

    [~, ~, score] = ...
        matvFRieszlike(Sigma_, n, nu, X);

    avg_n = mean(n);
    avg_nu = mean(nu);
    p = size(Omega_,1);
    c_1 = (avg_n^2*(avg_nu-p-2) + 2*avg_n)/((avg_nu-p)*(avg_nu-p-1)*(avg_nu-p-3));
    c_2 = (avg_n*(avg_nu-p-2)+avg_n^2+avg_n)/((avg_nu-p)*(avg_nu-p-1)*(avg_nu-p-3));
    c_4 = (avg_n-p-1)/((avg_n+avg_nu-1)*(avg_n+avg_nu+2))*((avg_n-p-2+1/(avg_nu+avg_n))*c_2-(1+(avg_n-p-1)/(avg_n+avg_nu))*c_1);
    c_3 = (avg_n-p-1)/(avg_n+avg_nu)*((avg_n-p-2)*c_2 - c_1)-(avg_n+avg_nu+1)*c_4;
    ckron2 = (avg_nu-(avg_n+avg_nu)*(c_3+c_4));
    cvec2 = (avg_n+avg_nu)*c_4;

    G = Dmatrix(p);
    invSig = inv(Sigma_);
    fisherinfo_Sigma_F = 1/2*G'*(ckron2*kron2(invSig) - cvec2*vec2(invSig))*G;
    
    for ii = 1:size(X,3)
        score.rc_paper(:,:,ii) = ...
            ivech(mean(diag(Y))*(fisherinfo_Sigma_F\score.Sigma_(ii,:)'));
    end

    varargout{1} = score;

end

if nargout >= 5

    [~, ~, ~, hessian, param] = ...
        matvFRieszlike(Sigma_, n, nu, X);
    
    varargout{2} = hessian;
    varargout{3} = param;
end

%%
% p = size(Omega_,1);
% 
% [~,iG] = Dmatrix(p);
% I = speye(p);
% L = ELmatrix(p);
% 
% Y = matvFRieszexpmat(n,nu);
% C_Om = chol(Omega_,'lower');
% C_Sig = C_Om/sqrtm(Y);
% Sigma_ = C_Sig*C_Sig';
% dOmega_dSigma = iG*kron(C_Sig*Y,I)*L'/(iG*kron(C_Sig,I)*L');
% 
% [nLogL, logLcontr, score, ~, ~, fisherinfo] = ...
%     matvFRieszlike(Sigma_, n, nu, X);
% 
% score.Omega_scaledbyiFish = ...
%     ivech(dOmega_dSigma*(fisherinfo.Sigma_\score.Sigma_'));

end
