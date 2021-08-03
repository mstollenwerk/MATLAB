function A = matvFrnd( Sigma_, df_1, df_2, N )
%MATVFRND
%
% REFERENCE: 
%   [1] Mulder & Pericci (2017) - 
%   [2] Gupta & Nagar (2000) - Theorem 5.2.2, Theorem 5.2.5, Theorem 5.3.2

p = size(Sigma_,1);
A = NaN(p,p,N);

sqrtSigma_ = sqrtm(Sigma_);
for ii=1:N

    S_1 = wishrnd( eye(p), df_1 );
    S_2 = wishrnd( eye(p), df_2 + p - 1 );
    
    sqrtinvS_2 = sqrtm(inv(S_2));
    
    A_ = sqrtSigma_*sqrtinvS_2*S_1*sqrtinvS_2*sqrtSigma_;
    
    A(:,:,ii) = sym_covmat(A_); % To make it exactly symmetric again (which is not given anymore after sqrtm.
        
end

end

