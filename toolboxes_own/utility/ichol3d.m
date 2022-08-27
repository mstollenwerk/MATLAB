function [Y] = ichol3d(X)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
Y = NaN(size(X));
for ii = 1:size(X,3)
    Y(:,:,ii) = X(:,:,ii)*X(:,:,ii)';
end
end

