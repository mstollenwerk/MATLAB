function R = matvsRieszrnd( Sigma_, n, T)
%MATVSRIESZRND
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
p = size(Sigma_,1);
C = chol(Sigma_,'lower');
M = diag(sqrt(1./n));

R = NaN(p,p,T);
for tt = 1:T
    
    BL = BarlettL(n);
    L = C*M*BL;
    R(:,:,tt) = L*L';
    
end

end

