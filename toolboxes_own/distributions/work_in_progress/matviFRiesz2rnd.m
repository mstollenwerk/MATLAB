function R = matviFRiesz2rnd( Omega_, n, nu, T)
%
p = length(n);
C = chol(Omega_,'lower');

R = NaN(p,p,T);
for ii = 1:T
    
    BL = BarlettL(n);
    BU = BarlettU(nu);
    
    L = C*BL/BU';
    
    R(:,:,ii) = L*L';
    
end

end

