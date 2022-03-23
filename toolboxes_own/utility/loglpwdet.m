function [y] = loglpwdet(X,powers,varargin)
%LPWDET Lower Power Weighted Determinant
%   Detailed explanation goes here
if size(powers,1)<size(powers,2)
    powers = powers';
end
if nargin >= 3 % Allow direct input of diagL
    if size(varargin{1}) ~= size(powers)
        error('Third input must be column vector of diagonal entries of Cholesky decomposition')
    end
    diagL = varargin{1};
else
    diagL = diag(chol(X,'lower'));
end

y = sum((2*powers).*log(diagL));

% powers(end+1) = 0;
% 
% y = 1;
% for ii = 1:p
%     y = y*det(X(1:ii,1:ii))^(powers(ii)-powers(ii+1));
% end

end

