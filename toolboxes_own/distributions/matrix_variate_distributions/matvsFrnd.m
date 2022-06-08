function R = matvsFrnd( Sigma_, n, nu, N )
%MATVSFRND Random standardized matrix-F distribution matrix.
%   Detailed explanation goes here

p = size(Sigma_,1);
Omega_ = (nu-p-1)/n*Sigma_;
R = matvFrnd( Omega_, n, nu, N );

end

