function A = matviFRiesz2rnd( Sigma_, n, nu, T)
%
p = length(n);
X = matvRieszrnd(eye(p), nu, T);
Y = matviRiesz2rnd(eye(p), n, T);
C = chol(Sigma_,'lower');

A = NaN(p,p,T);
for ii = 1:T
    Cx = chol(X(:,:,ii),'lower');
    A(:,:,ii) = C*Cx*Y(:,:,ii)*Cx'*C';
end

end

