function y = tril3d(A,k)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
[p,~,N] = size(A);
id_ = tril(true(p),k);

y = NaN(N,sum(id_(:)));
for ii = 1:N
    Xii = A(:,:,ii);
    y(ii,:) = Xii(id_);
end

end