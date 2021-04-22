function y = mvgammaln(x,d)
%MVGAMMALC computes the natural logarithm of the multivariate gamma function 
%
% USAGE:
%  Y = mvgammaln(x,d)
%
% INPUTS:
%   X            - x-value
%   D            - "number of sums" parameter
%
% OUTPUTS:
%   Y            - y-value
%
% COMMENTS:
%   Used in the probability density function of the Wishart and inverse 
%   Wishart distributions.
%   Gamma_d(x) = pi^(d(d-1)/4) \prod_(j=1)^d Gamma(x+(1-j)/2)
%   log(Gamma_d(x)) = d(d-1)/4 log(pi) + \sum_(j=1)^d log(Gamma(x+(1-j)/2))
%
% REFERENCES:
%      [1] James (1964) - Distributions of Matrix Variates and Latent 
%      Roots Derived from Normal Samples. 

% Michael Stollenwerk
% michael.stollenwerk@live.com
% 12.02.2017

y = d*(d-1)/4*log(pi)+sum(gammaln(x+(1-(1:d))/2));

end