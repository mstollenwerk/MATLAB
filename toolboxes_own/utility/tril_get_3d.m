function [Out_] = tril_get_3d(Matrix,RowAboveMainDiag)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
[p,~,N] = size(Matrix);

sel = tril(true(p),RowAboveMainDiag);
dim_out = sum(sum(sel));

Out_ = NaN(N,dim_out);
for ii = 1:N
    Mat_ii = Matrix(:,:,ii);
    Out_(ii,:) = Mat_ii(sel);
end

end