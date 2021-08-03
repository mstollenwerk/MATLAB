function A = matvRiesz2rnd( Sigma_, n, T)
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
p_ = p*(p+1)/2 - p;
U = chol(Sigma_);

A = NaN(p,p,T);
for tt = 1:T
    chi2_ = NaN(p,1);
    for pp = 1:p
        chi2_(pp) = sqrt(chi2rnd( n(pp) - p + pp ));
    end
    H = diag(chi2_);
    H(triu(true(p),1)) = randn(p_,1);
    A_ = U*(H*H')*U';
    A(:,:,tt) = (A_ + A_')./2; % ensure symmetry
end

end

