function [ L, d ] = dpceig( X, L_ )
%DPCEIG sorts eigenvalues big2small (eigenvectors accordingly). Diagonal
%elements of eigenvectors matrix are restricted to be positive (1 input) or
%are restricted s.th. diagonal elements of L_*X are positive (2 inputs).

% Michael Stollenwerk
% michael.stollenwerk@live.com
% 29.08.2016

[L, d]=eig(X,'vector');
[d,perm]=sort(d,'descend');
L=L(:,perm);
if nargin==1
    L = bsxfun(@times,L,sign(diag(L)'));
elseif nargin==2
    L = bsxfun(@times,L,sign(diag(L_*L)'));
end

end