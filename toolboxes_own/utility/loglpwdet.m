function [y] = loglpwdet(X,powers)
%LPWDET Lower Power Weighted Determinant
%   Detailed explanation goes here
if size(powers,1)<size(powers,2)
    powers = powers';
end
diagL = diag(chol(X,'lower'));

y = sum((2*powers).*log(diagL));

% powers(end+1) = 0;
% 
% y = 1;
% for ii = 1:p
%     y = y*det(X(1:ii,1:ii))^(powers(ii)-powers(ii+1));
% end

end

