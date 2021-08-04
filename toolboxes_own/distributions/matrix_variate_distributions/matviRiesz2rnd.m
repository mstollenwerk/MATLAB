function A = matviRiesz2rnd( Sigma_, n, T)
%MATVRIESZ2RND
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

invA = matvRiesz2rnd(inv(Sigma_), n, T);

A = NaN(p,p,T);
for tt = 1:T
    A(:,:,tt) = inv(invA(:,:,tt));
end

end

