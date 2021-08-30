function [y] = logupwdet(X,powers,varargin)
%LPWDET Lower Power Weighted Determinant
%   Detailed explanation goes here

% REFERENCES:
%      [1] Blasques, Lucas, Opschoor and Rossini (2021)

% Michael Stollenwerk
% michael.stollenwerk@live.com
% 30.08.2021

if size(powers,1)<size(powers,2)
    powers = powers';
end
if nargin >= 3 % Allow direct input of diagL
    diagU = varargin{1};
else
    diagU = diag(inv(chol(inv(X), 'lower')'));
end

y = sum((2*powers).*log(diagU));

end