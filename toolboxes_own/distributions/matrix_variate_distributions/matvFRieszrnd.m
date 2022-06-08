function R = matvFRieszrnd( Omega_, n, nu, N)
%
p = length(n);
C = chol(Omega_,'lower');

R = NaN(p,p,N);
for ii = 1:N
    
    BL = BarlettL(n);
    BU = BarlettU(nu);
    
    L = C/BU'*BL;
    
    R(:,:,ii) = L*L';
    
end

end

