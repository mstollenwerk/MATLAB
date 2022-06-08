function R = matvsitRiesz2rnd( Sigma_, n, nu, N )
%MATVSITRIESZ2RND Random standardized inverse t-Riesz 2 matrix.
%

p = size(Sigma_,1);
C = chol(Sigma_,'lower');
m = diag(matviRieszexpmat(nu));
M = diag(sqrt(1./m));

R = NaN(p,p,N);
for ii = 1:N
    
    BU = BarlettU(nu);
    L = C*M/BU';
    siR2 = L*L';
    
    gam = gamrnd( n/2, 2/n );
    
    R(:,:,ii) = siR2.*gam;
    
%     chi2_ = chi2rnd(nu); 
%     
%     R(:,:,ii) = Riesz_/(chi2_/nu);

end

end
