function [ nLogL, logLcontr, varargout ] = matvsWishlike( Omega_, df, X, varargin )
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

if nargout <= 2
    
    [nLogL, logLcontr] = ...
        matvWishlike(Sigma_, df, X);
    
elseif nargout >= 3

    [nLogL, logLcontr, score, ~, ~, fisherinfo] = ...
        matvWishlike(Sigma_, df, X);

    score.Omega_scaledbyiFish = ...
        ivech(dOmega_dSigma*(fisherinfo.Sigma_\score.Sigma_'));

    score.rc_paper = score.Omega_scaledbyiFish;

    varargout{1} = score;
end

end
