function y = lgmvgammaln(x)
%MVGAMMALC computes the natural logarithm of the lower generalized multivariate gamma function 
%
% USAGE:
%  Y = lgmvgammaln(x,d)
%
% INPUTS:
%   X            - x-value
%   D            - "number of sums" parameter
%
% OUTPUTS:
%   Y            - y-value
%
% COMMENTS:
%
% REFERENCES:
%      [1] Blasques, Lucas, Opschoor and Rossini (2021)

% Michael Stollenwerk
% michael.stollenwerk@live.com
% 11.02.2021

if any(size(x)==1)
    if size(x,1) > size(x,2)
        x = x';
    end
end
d = size(x,2);
y = d*(d-1)/4*log(pi)+sum(gammaln(x+(1-(1:d))/2),2);

end