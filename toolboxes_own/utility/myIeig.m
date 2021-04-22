function [ V, D ] = myIeig( X, I )
%Return eigenvector matrix with elements in indice vector "I" being
%positive
[V,D] = eig(X);
V = bsxfun(@times,V,sign(V(I)));
end

