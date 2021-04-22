function y = dcholX_dX(X)
%Tested. See
%https://mathoverflow.net/questions/150427/the-derivative-of-the-cholesky-factor
%three consecutive comments by pete in answer by Steven Pav.

p = size(X,1);

L = chol(X,'lower');
invL = inv(L);

EL = sparse(Ematrix(p)); 
D = sparse(Dmatrix(p));

Q = .5*(D*EL)'*D*EL;

y = EL*kron(eye(p),L)*Q*kron(invL,invL)*D;

end