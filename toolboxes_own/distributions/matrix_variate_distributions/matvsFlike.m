function [ nLogL, logLcontr, varargout ] = matvsFlike( Omega_, n, nu, X, varargin )
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

Y = n*nu/(nu-k-1);
Sigma_ = Omega_/Y;
dOmega_dSigma = Y;

[nLogL, logLcontr] = ...
    matvFlike(Sigma_, n, nu, X);
    
if nargout >= 3
    
    [nLogL, logLcontr, score, ~, param, fisherinfo] = ...
        matvFlike(Sigma_, n, nu, X);

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
