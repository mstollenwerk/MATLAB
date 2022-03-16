function [ nLogL, logLcontr, score ] = matvsiRieszlike( Omega_, n, X, varargin )
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

Y = matviRieszexpmat(n);
iCdot_Om = inv(chol(inv(Omega_),'lower'));
iCdot_Sig = iCdot_Om/sqrtm(Y);
Sigma_ = iCdot_Sig'*iCdot_Sig;
dOmega_dSigma = iG*kron(iCdot_Sig',iCdot_Sig'*Y*iCdot_Sig)*L'/(iG*kron(iCdot_Sig',Sigma_)*L');

if nargout <= 2
    
    [nLogL, logLcontr] = ...
        matviRieszlike(Sigma_, n, X);
    
elseif nargout >= 3

    [nLogL, logLcontr, score, ~, ~, fisherinfo] = ...
        matviRieszlike(Sigma_, n, X);

    score.Omega_scaledbyiFish = ...
        ivech(dOmega_dSigma*(fisherinfo.Sigma_\score.Sigma_'));

    varargout{1} = score;

end

end
