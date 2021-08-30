function [y] = upwdet(X,powers)
%UPWDET Upper Power Weighted Determinant
%   Detailed explanation goes here
if size(powers,1)<size(powers,2)
    powers = powers';
end


U = inv(chol(inv(X), 'lower')');
diagU = diag(U);
y = prod(diagU.^(2*powers));
end

