function A = matvsiRiesz2rnd( Sigma_, n, T)
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
p = size(Sigma_,1);
p_ = p*(p+1)/2 - p;
C = chol(Sigma_,'lower');

A = NaN(p,p,T);
for tt = 1:T
    chi2_ = NaN(p,1);
    for pp = 1:p
        chi2_(pp) = sqrt(chi2rnd( n(pp) + pp - p ));
    end
    G = diag(chi2_);
    G(triu(true(p),1)) = randn(p_,1);
    iG = inv(G);
    A_ = C*(iG'*iG)*C';
    A(:,:,tt) = (A_ + A_')./2; % ensure symmetry
end

end

