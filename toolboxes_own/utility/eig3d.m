function e = eig3d(X)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
[p,~,N] = size(X);

e = NaN(N,p);
for ii = 1:N
    e(ii,:) = eig(X(:,:,ii));
end

end

