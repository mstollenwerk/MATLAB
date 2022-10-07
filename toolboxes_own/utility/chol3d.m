function [Y] = chol3d(X)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
Y = NaN(size(X));
for ii = 1:size(X,3)
    Y(:,:,ii) = chol(X(:,:,ii),'lower');
end
end

