function logdetX = logdet(X)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

[L, flag] = chol(X);
if flag == 0 
    logdetX = 2*sum(log(diag(L)));
else
    logdetX = log(det(X));
end

end

