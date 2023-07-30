function R = matvsiFRiesz2rnd( Sigma_, n, nu, N)
%MATVSRND Random standardized inverse F-Riesz 2 matrices.
%   Detailed explanation goes here
%
% Michael Stollenwerk
% michael.stollenwerk@live.com
% 06.06.2022
%
% DEPENDENCIES:
%
p = size(Sigma_,1);
C = chol(Sigma_,'lower');
m = diag(matviFRiesz2expmat(n,nu));
M = diag(sqrt(1./m));


R = NaN(p,p,N);
for ii = 1:N
    
    BL = BarlettL(n);    
    BU = BarlettU(nu);
    
    L = C*M*BL/BU';
    
    R(:,:,ii) = L*L';
    
end

end

