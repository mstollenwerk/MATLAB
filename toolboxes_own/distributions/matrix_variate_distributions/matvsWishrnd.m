function R = matvsWishrnd( Sigma_, df, T )
%MATVWISHRND Random Wishart matrix.
%
% Michael Stollenwerk
% michael.stollenwerk@live.com
% 14.08.2020

Omega_ = Sigma_./df;
R = matvWishrnd( Omega_, df, T );

end

