function y = diag3d(X)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
[p,~,N] = size(X);

y = NaN(N,p);
for ii = 1:N
    y(ii,:) = diag(X(:,:,ii));
end

end

