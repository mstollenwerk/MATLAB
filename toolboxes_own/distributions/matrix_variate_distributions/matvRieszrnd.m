function R = matvRieszrnd( Omega_, n, T)
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
% 11.02.2021
%
% DEPENDENCIES:
%
p = size(Omega_,1);
C = chol(Omega_,'lower');

R = NaN(p,p,T);
for tt = 1:T
    
    BL = BarlettL(n);
    L = C*BL;
    R(:,:,tt) = L*L';
    
end

end

