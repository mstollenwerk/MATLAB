function outputArg1 = matvcov(X)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
[p,~,n] = size(X);
Y = NaN(p^2,p^2,n);
meanX = mean(X,3);
for ii = 1:n
    vecMeanDeviation = vec(X(:,:,ii)-meanX);
    Y(:,:,ii) = vecMeanDeviation*vecMeanDeviation';
end

outputArg1 = mean(Y,3);
end

