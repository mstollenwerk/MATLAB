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
    
[nLogL, logLcontr] = ...
    matviWishlike(Sigma_, n, X);
    
if nargout >= 3

    [~, ~, score, ~, param, fisherinfo] = ...
        matviWishlike(Sigma_, n, X);
    
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
