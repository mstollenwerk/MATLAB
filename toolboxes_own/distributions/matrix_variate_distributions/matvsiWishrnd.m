function R = matvsiWishrnd( Sigma_, nu, T )
%MATVIWISHRND Random inverse-Wishart matrix.
%
% Michael Stollenwerk
% michael.stollenwerk@live.com

p = size(Sigma_,1);
Omega_ = (nu-p-1)*Sigma_;
R = matviWishrnd( Omega_, nu, T );

end