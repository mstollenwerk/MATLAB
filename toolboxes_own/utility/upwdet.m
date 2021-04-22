function [y] = upwdet(X,powers)
%UPWDET Upper Power Weighted Determinant
%   Detailed explanation goes here
diagL = diag(chol(X));
y = prod(diagL.^(2*powers));
end

