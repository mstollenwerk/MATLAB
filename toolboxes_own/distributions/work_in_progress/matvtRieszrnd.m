function R = matvtRieszrnd( Omega_, n, nu, N )
%MATVTRIESZRND Random t-Riesz 1 matrix
%

p = size(Omega_,1);
C = chol(Omega_,'lower');

R = NaN(p,p,N);
for ii = 1:N
    
    BL = BarlettL(n);
    L = C*BL;
    Riesz_ = L*L';
    
    chi2_ = chi2rnd(nu); 
    
    R(:,:,ii) = Riesz_/(chi2_/nu);
    
    % Alternatively:
%     gamma_ = gamrnd( nu/2, 2/nu );
%     R(:,:,tt) = Riesz_./gamma_;

end

end

