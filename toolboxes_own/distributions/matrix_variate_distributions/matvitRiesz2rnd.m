function R = matvitRiesz2rnd( Omega_, n, nu, N )
%MATVITRIEZS2RND Random Inverse t-Riesz 2 matrix.
%   Detailed explanation goes here

p = size(Omega_,1);
C = chol(Omega_,'lower');

R = NaN(p,p,N);
for ii = 1:N
    
    BU = BarlettU(nu);
    L = C/BU';
    iR2 = L*L';
    
    chi2_ = chi2rnd(n); 
    R(:,:,ii) = iR2*(chi2_/n);

    % Alternatively:
%     gam = gamrnd( n/2, 2/n );
%     R(:,:,ii) = iR2.*gam;
    
end

end

