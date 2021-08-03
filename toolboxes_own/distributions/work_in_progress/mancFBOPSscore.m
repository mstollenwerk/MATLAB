function [ score ] = mancFBOPSscore( X, Sigma_, Omega_, nu1, nu2 )
%MANCFBOPSSCORE Score for the rescaled matrix-variate non-central F 
%   as defined by Bodnar, Okhrin, Parolya and Stollenwerk (2020).
%
% USAGE:
% 
%
% INPUTS:
%   SIGMA_ - (p by p) regulates the covariance matrix.
%   OMEGA_ - (p by p) non-centrality matrix.
%   NU1    - Double, first degrees of freedom parameter.
%   NU2    - Double, second degrees of freedom parameter.
%
% OUTPUTS:
%   SCORE        - Struct with fields as param names. Contains derivatives
%                  of pdf w.r.t. parameters.
%
% COMMENTS:
%
% REFERENCES:
%      [1] Bodnar, Okhrin, Parolya and Stollenwerk (2020)
%
% DEPENDENCIES:
%
%
%  See also
%
% Michael Stollenwerk
% michael.stollenwerk@live.com
% 16.04.2020

p = size(X,1);

c_nu = (nu1+nu2)/2/(n2-p-1);
cc_nu = nu1/(n2-p-1);

Z = ...
    - nu1/2/Sigma_ + c_nu\(Sigma_ * (inv(X) + cc_nu\Sigma_) * Sigma_) ... % Opschoor part
    + 1/2*Sigma_\Omega_/Sigma_ - c_nu*(2*X\Sigma_*Omega_ + cc_nu*Omega_); % Non-central part

score.Sigma_ = 2*Z - Z.*eye(p);

end

