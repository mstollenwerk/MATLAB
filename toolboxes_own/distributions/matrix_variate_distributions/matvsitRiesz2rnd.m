function R = matvsitRiesz2rnd( Sigma_, n, nu, N )
%MATVSITRIESZ2RND Random standardized inverse t-Riesz 2 matrix.
%

Omega_ = matvStandardize('itRiesz2',Sigma_,[n; nu]);
R = matvitRiesz2rnd( Omega_, n, nu, N );

end
