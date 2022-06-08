function [ nLogL, logLcontr, score ] = matvsRiesz2like( Omega_, n, X, varargin )
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
L = ELmatrix(p);

Y = diag(n);
U_Om = cholU(Omega_);
U = U_Om/sqrtm(Y);
Sigma_ = U*U';
dOmega_dSigma = iG*kron(U,Omega_)*L'/(iG*kron(U,Sigma_)*L');

[nLogL, logLcontr, score, ~, ~, fisherinfo] = ...
    matvRiesz2like(Sigma_, n, X);

score.Omega_scaledbyiFish = ...
    ivech(dOmega_dSigma*(fisherinfo.Sigma_\score.Sigma_'));

end
