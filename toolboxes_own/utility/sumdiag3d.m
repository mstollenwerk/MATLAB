function y = sumdiag3d(X)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
y = NaN(size(X,3),1);
for ii = 1:size(X,3)
    y(ii) = sum(diag(X(:,:,ii)));
end
end

