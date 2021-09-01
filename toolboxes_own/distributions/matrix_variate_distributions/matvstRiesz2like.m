function [ nLogL, logLcontr, score ] = matvstRiesz2like( Omega_, n, nu, X, varargin )
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
% 31.08.2021

p = size(Omega_,1);

[~,iG] = Dmatrix(p);
L = ELmatrix(p);

Y = nu/(nu-2)*diag(n);
U_Om = chol(Omega_);
U = U_Om/sqrtm(Y);
Sigma_ = U*U';
dOmega_dSigma = iG*kron(U,Omega_)*L'/(iG*kron(U,Sigma_)*L');

[nLogL, logLcontr, score, ~, ~, fisherinfo] = ...
    matvtRiesz2like(Sigma_, n, nu, X);

score.Omega_scaledbyiFish = ...
    ivech(dOmega_dSigma*(fisherinfo.Sigma_\score.Sigma_'));

end
