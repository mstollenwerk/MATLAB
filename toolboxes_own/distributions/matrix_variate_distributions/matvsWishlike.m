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

[nLogL, logLcontr] = ...
    matvWishlike(Sigma_, df, X);
    
if nargout >= 3

    [nLogL, logLcontr, score, ~, param, fisherinfo] = ...
        matvWishlike(Sigma_, df, X);

    for ii = 1:size(X,3)
        score.Omega_scaledbyiFish(:,:,ii) = ...
            ivech(dOmega_dSigma*(fisherinfo.Sigma_\score.Sigma_(ii,:)'));
    end

    score.rc_paper = score.Omega_scaledbyiFish;

    varargout{1} = score;
    varargout{3} = param;
    varargout{4} = fisherinfo;

end

end
