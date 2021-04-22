function y = mvgamma(x,d)
%MVGAMMALC computes the multivariate gamma function 
%
% USAGE:
%  Y = mvgamma(x,d)
%
% INPUTS:
%   X            - x-value
%   D            - "number of products" parameter
%
% OUTPUTS:
%   Y            - y-value
%
% COMMENTS:
%   Used in probability density functions  (e.g. (Inverse) Wishart)
%   Gamma_d(x) = pi^(d(d-1)/4) \prod_(j=1)^d Gamma(x+(1-j)/2)
%
% REFERENCES:
%      [1] James (1964) - Distributions of Matrix Variates and Latent 
%      Roots Derived from Normal Samples. 

% Michael Stollenwerk
% michael.stollenwerk@live.com
% 17.04.2019

y = pi^(d*(d-1)/4) * prod(gamma(x+(1-(1:d))/2));

end
