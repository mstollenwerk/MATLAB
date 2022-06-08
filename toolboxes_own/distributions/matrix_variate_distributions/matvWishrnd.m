function R = matvWishrnd( Omega_, df, T )
%MATVWISHRND Random Wishart matrix.
%
% Michael Stollenwerk
% michael.stollenwerk@live.com
% 14.08.2020

p = size(Omega_,1);
C = chol(Omega_,'lower');

R = NaN(p,p,T);
for ii=1:T
    
    B = BarlettL(ones(p,1)*df);
    L = C*B;
    R(:,:,ii) = L*L';
    
end

end

