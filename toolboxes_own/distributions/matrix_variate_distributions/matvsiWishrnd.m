function R = matvsiWishrnd( Sigma_, df, T )
%MATVIWISHRND Random inverse-Wishart matrix.
%
% Michael Stollenwerk
% michael.stollenwerk@live.com

p = length(df);
Omega_ = (n-p-1)*Sigma_;
R = matviWishrnd( Omega_, df, T );

end