function logdetX = logdet3d(X)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
logdetX = NaN(size(X,3),1);
for ii = 1:size(X,3)
    if all(isnan(squeeze(X(:,:,ii))))
        logdetX(ii) = NaN;
    else
        logdetX(ii) = logdet(X(:,:,ii));
    end
end

end