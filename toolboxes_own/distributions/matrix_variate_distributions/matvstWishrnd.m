function R = matvstWishrnd( Sigma_, n, nu, T )

Omega_ = Sigma_*(nu-2)/nu/n;
R = matvtWishrnd( Omega_, n, nu, T );

end
