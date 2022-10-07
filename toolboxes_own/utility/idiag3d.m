function y = idiag3d(X)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

[N,p] = size(X);

y = NaN(p,p,N);
for ii = 1:N
    y(:,:,ii) = diag(X(ii,:));
end 

end
