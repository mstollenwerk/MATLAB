function E = Ematrix(n)
%EMATRIX create the elimination matrix as in Harville - Matrix Algebra from
%a statisticians perspective p. 358 with c_i = 1;
%
% Michael Stollenwerk
% michael.stollenwerk@live.com
% 25.06.2018

E = eye(n^2);
idx_delete_row = tril(true(n));
E = E(idx_delete_row(:),:);

end