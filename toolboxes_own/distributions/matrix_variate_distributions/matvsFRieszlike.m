function [ nLogL, logLcontr, score ] = matvsFRieszlike( Omega_, n, nu, X, varargin )
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

p = size(Omega_,1);

[~,iG] = Dmatrix(p);
I = speye(p);
L = ELmatrix(p);

Y = matvFRieszexpmat(n,nu);
C_Om = chol(Omega_,'lower');
C_Sig = C_Om/sqrtm(Y);
Sigma_ = C_Sig*C_Sig';
dOmega_dSigma = iG*kron(C_Sig*Y,I)*L'/(iG*kron(C_Sig,I)*L');

[nLogL, logLcontr, score, ~, ~, fisherinfo] = ...
    matvFRieszlike(Sigma_, n, nu, X);

score.Omega_scaledbyiFish = ...
    ivech(dOmega_dSigma*(fisherinfo.Sigma_\score.Sigma_'));

end
