function R = matvsiRiesz2rnd( Sigma_, nu, T)
%MATVSIRIESZ2RND Random standardized inverse Riesz 2 matrix.
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
% 06.06.2022
%
% DEPENDENCIES:
%

p = size(Sigma_,1);
C = chol(Sigma_,'lower');
m = diag(matviRiesz2expmat(nu));
M = diag(sqrt(1./m));

R = NaN(p,p,T);
for tt = 1:T
    
    BU = BarlettU(nu);
    L = C*M/BU';
    R(:,:,tt) = L*L';
    
end

end

