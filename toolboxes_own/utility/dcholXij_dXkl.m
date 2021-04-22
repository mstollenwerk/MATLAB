function y = dcholXij_dXkl(X,i,j,k,l)
%Tested. See
%https://mathoverflow.net/questions/150427/the-derivative-of-the-cholesky-factor
%answer by Ian Murray.

p = size(X,1);

L = chol(X,'lower');
Li = inv(L);

term1 = 0;
for m = (j+1):p
    term1 = term1 + L(i,m)*Li(m,k);
end
term2 = .5*L(i,j)*Li(j,k);

y = (term1 + term2)*Li(j,l);

if k ~= l
    term3 = 0;
    for m = (j+1):p
        term3 = term3 + L(i,m)*Li(m,l);
    end
    term4 = .5*L(i,j)*Li(j,l);
    
    y = y + (term3 + term4)*Li(j,k);
end

end