function R = matviRiesz2rnd( Omega_, nu, T)
%MATVRIESZRND
%
% USAGE:
%   
%
% INPUTS:
%   
%
% OUTPUTS:
%   
%  See also 
%
% COMMENTS:
% 
% REFERENCES:
%      [1] Stollenwerk (2020)
%
% Michael Stollenwerk
% michael.stollenwerk@live.com
% 24.03.2022
%
% DEPENDENCIES:
%

p = size(Omega_,1);
C = chol(Omega_,'lower');

R = NaN(p,p,T);
for tt = 1:T
    
    BU = BarlettU(nu);
    L = C/BU';
    R(:,:,ii) = L*L';
    
end

end

