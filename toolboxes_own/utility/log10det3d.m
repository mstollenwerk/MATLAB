function log10detX = log10det3d(X)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
log10detX = NaN(size(X,3),1);
for ii = 1:size(X,3)
    L = chol(X(:,:,ii));
    log10detX(ii) = 2*sum(log10(diag(L)));
end

end