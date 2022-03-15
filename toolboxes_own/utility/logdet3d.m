function logdetX = logdet3d(X)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
logdetX = NaN(size(X,3),1);
for ii = 1:size(X,3)
    if all(isnan(squeeze(X(:,:,ii))))
        logdetX(ii) = NaN;
    else
        L = chol(X(:,:,ii));
        logdetX(ii) = 2*sum(log(diag(L)));
    end
end

end