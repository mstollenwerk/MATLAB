function [ nLogL ] = matvsFRiesz2like( Omega_, n, nu, X, varargin )
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
% 08.09.2021

Y = matvFRiesz2expmat(n,nu);
C_Om = cholU(Omega_);
C_Sig = C_Om/sqrtm(Y);
Sigma_ = C_Sig*C_Sig';
nLogL = matvFRieszlike(Sigma_, n, nu, X);

% p = size(Omega_,1);
% [nLogL, logLcontr, score, ~, ~, fisherinfo] = ...
%     matvFRieszlike(Sigma_, n, nu, X);
% [~,iG] = Dmatrix(p);
% I = speye(p);
% L = ELmatrix(p);
% dOmega_dSigma = NaN;
% score.Omega_scaledbyiFish = ...
%     ivech(dOmega_dSigma*(fisherinfo.Sigma_\score.Sigma_'));

end
