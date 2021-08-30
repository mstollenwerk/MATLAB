function [ nLogL, logLcontr, score ] = matvstWishlike( Omega_, n, nu, X, varargin )
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

Y = n*nu/(nu-2);
Sigma_ = Omega_/Y;
dOmega_dSigma = Y;

[nLogL, logLcontr, score, ~, ~, fisherinfo] = ...
    matvtWishlike(Sigma_, n, nu, X);

score.Omega_scaledbyiFish = ...
    ivech(dOmega_dSigma*(fisherinfo.Sigma_\score.Sigma_'));

end