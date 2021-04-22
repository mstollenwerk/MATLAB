function invX = inv3d(X)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
invX = NaN(size(X));
for ii = 1:size(X,3)
    invX(:,:,ii) = inv(X(:,:,ii));
end
end

