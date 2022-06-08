function R = matvsiWishrnd( Sigma_, df, T )
%MATVIWISHRND Random inverse-Wishart matrix.
%
% Michael Stollenwerk
% michael.stollenwerk@live.com

p = size(Sigma_,1);
Omega_ = (n-p-1)*Sigma_;
R = matviWishrnd( Omega_, df, T );

end