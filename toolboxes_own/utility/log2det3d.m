function log2detX = log2det3d(X)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
log2detX = NaN(size(X,3),1);
for ii = 1:size(X,3)
    L = chol(X(:,:,ii));
    log2detX(ii) = 2*sum(log2(diag(L)));
end

end