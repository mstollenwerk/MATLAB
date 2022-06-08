function R = matvFrnd( Omega_, n, nu, N )
%MATVFRND Random matrix-F distribution matrix.
%   Detailed explanation goes here

p = length(n);
C = chol(Omega_,'lower');

R = NaN(p,p,N);
for ii = 1:N
    
    BL = BarlettL(ones(p,1).*n);
    BU = BarlettU(ones(p,1).*nu);
    
    L = C/BU'*BL;
    
    R(:,:,ii) = L*L';
    
end


end

