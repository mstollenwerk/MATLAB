function [Y] = det_(X)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
N = size(X,3);

Y = NaN(N,1);
for ii = 1:N
    Y(ii) = det(X(:,:,ii));
end

end

