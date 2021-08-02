function trX = trace3d(X)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
trX = NaN(size(X,3),1);
for ii = 1:size(X,3)
    trX(ii) = trace(X(:,:,ii));
end
end

