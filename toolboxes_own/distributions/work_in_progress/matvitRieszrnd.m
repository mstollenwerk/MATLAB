function A = matvitRieszrnd( Sigma_, n, nu, N )
%MATVTWISHRND
%

A = matviRieszrnd( Sigma_, n, N );
for ii=1:N

    chi2_ = chi2rnd(nu); 
    A(:,:,ii) = A(:,:,ii).*(chi2_/nu);
    
    % Alternatively:
%     gamma_ = gamrnd( df_t/2, 2/df_t );
%     A(:,:,ii) = A(:,:,ii).*gamma_;
    
end

end

