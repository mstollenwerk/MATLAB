function [ nLogL, logLcontr, score ] = matvsiWishlike( Omega_, df, X, varargin )
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

Y = 1/(n-k-1);
Sigma_ = Omega_/Y;
dOmega_dSigma = Y;

[nLogL, logLcontr, score, ~, ~, fisherinfo] = ...
    matviWishlike(Sigma_, df, X);

score.Omega_scaledbyiFish = ...
    ivech(dOmega_dSigma*(fisherinfo.Sigma_\score.Sigma_'));

end
