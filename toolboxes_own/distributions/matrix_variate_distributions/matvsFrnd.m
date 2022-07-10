function R = matvsFrnd( Sigma_, n, nu, N )
%MATVSFRND Random standardized matrix-F distribution matrix.
%   Detailed explanation goes here

Omega_ = matvStandardize('F',Sigma_,[n,nu]);
R = matvFrnd( Omega_, n, nu, N );

end

