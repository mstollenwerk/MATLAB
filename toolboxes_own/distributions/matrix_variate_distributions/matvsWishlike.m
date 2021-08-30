function [ nLogL, logLcontr, score ] = matvsWishlike( Omega_, df, X, varargin )
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

Sigma_ = Omega_/df;
dOmega_dSigma = df;

[nLogL, logLcontr, score, ~, ~, fisherinfo] = ...
    matvWishlike(Sigma_, df, X);

score.Omega_scaledbyiFish = ...
    ivech(dOmega_dSigma*(fisherinfo.Sigma_\score.Sigma_'));

end
