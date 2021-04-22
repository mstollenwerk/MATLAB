function z = mvbetaln(x,y,d)
%MVGAMMALC Natural logarithm of multivariate Beta function.
%
% USAGE:
%  Y = mvbetaln(x,y,d)
%
% INPUTS:
%   X            - Double
%   Y            - Double
%   D            - Double, "number of sums" parameter
%
% OUTPUTS:
%   Z            - Double
%
% COMMENTS:
%   Used in the probability density function of the Wishart and inverse 
%   Wishart distributions.
%   Gamma_d(x) = pi^(d(d-1)/4) \prod_(j=1)^d Gamma(x+(1-j)/2)
%   log(Gamma_d(x)) = d(d-1)/4 log(pi) + \sum_(j=1)^d log(Gamma(x+(1-j)/2))
%   Beta_d(x,y) = (Gamma_d(x) * Gamma_d(y))/Gamma_d(x+y)
%   log(Beta_d(x,y)) = log(Gamma_d(x)) + log(Gamma_d(y)) - log(Gamma_d(x+y))
%
% REFERENCES:
%      [1] https://dlmf.nist.gov/35.3#E7
%
% Michael Stollenwerk
% michael.stollenwerk@live.com
% 16.04.2020

z = mvgammaln(x,d) ...
    + mvgammaln(y,d) ...
    - mvgammaln(x+y,d);

end