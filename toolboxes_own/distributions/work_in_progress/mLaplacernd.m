function [ mLaplace_rnd ] = mLaplacernd( Sigma_, Phi_ )
%MALRND Generate vector-valued (symmetric) Laplace random sample
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
%      [1] Kotz, Kozubowski and Podgorski (2001) - The Laplace distribution
%          and generalizations.
%
% Michael Stollenwerk
% michael.stollenwerk@live.com
% 14.02.2020
%
% DEPENDENCIES:
%
%% 
[p, ~, k] = size(Sigma_);
[n, ~, k_] = size(Phi_);
if k~=k_
    error('Input array dimensions do not agree.')
end

mLaplace_rnd = NaN(p,n,k);
for ii = 1:k
    W = exprnd(1);
    Sigma_vecX = kron(Sigma_(:,:,ii), Phi_(:,:,ii));
    vecX = mvnrnd(zeros(p*n,1), Sigma_vecX);
    mLaplace_rnd(:,:,ii) = reshape(sqrt(W)*vecX,n,p)';
end

end
