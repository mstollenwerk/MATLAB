function R = matviWishrnd( Omega_, df, T )
%MATVIWISHRND Random inverse-Wishart matrix.
%
% Michael Stollenwerk
% michael.stollenwerk@live.com

p = size(Omega_,1);
C = chol(Omega_,'lower');

R = NaN(p,p,T);
for ii=1:T

    BU = BarlettU(ones(p,1)*df);
    L = C/BU';
    R(:,:,ii) = L*L';
    
end

end

