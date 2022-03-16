function [ nLogL, logLcontr, varargout ] = matvsiWishlike( Omega_, n, X, varargin )
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

k = size(Omega_,1);

Y = 1/(n-k-1);
Sigma_ = Omega_/Y;
dOmega_dSigma = Y;

if nargout <= 2
    
    [nLogL, logLcontr] = ...
        matviWishlike(Sigma_, n, X);
    
elseif nargout >= 3

    [nLogL, logLcontr, score, ~, ~, fisherinfo] = ...
        matviWishlike(Sigma_, n, X);

    score.Omega_scaledbyiFish = ...
        ivech(dOmega_dSigma*(fisherinfo.Sigma_\score.Sigma_'));

    score.rc_paper = score.Omega_scaledbyiFish;

    varargout{1} = score;
end

end
