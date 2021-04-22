function logdetX = logdet(X)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
L = chol(X);
logdetX = 2*sum(log(diag(L)));

end

