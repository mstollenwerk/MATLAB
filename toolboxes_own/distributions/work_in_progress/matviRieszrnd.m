function A = matviRieszrnd( Sigma_, n, T)
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
p = size(Sigma_,1);

invA = matvRieszrnd(inv(Sigma_), n, T);

A = NaN(p,p,T);
for tt = 1:T
    A(:,:,tt) = inv(invA(:,:,tt));
end

end

