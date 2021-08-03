function A = matvFRieszrnd( Sigma_, n, nu, T)
%
p = length(n);
X = matvRieszrnd(eye(p), n, T);
Y = matviRiesz2rnd(eye(p), nu, T);
L_Sig = chol(Sigma_,'lower');

A = NaN(p,p,T);
for ii = 1:T
    L_Y = chol(Y(:,:,ii),'lower');
    A(:,:,ii) = L_Sig*L_Y*X(:,:,ii)*L_Y'*L_Sig';
end

end

