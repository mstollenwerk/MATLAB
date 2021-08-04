function [ vLaplace_rnd ] = vLaplacernd( Sigma_ )
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
[p, ~, K] = size(Sigma_);

vLaplace_rnd = NaN(p,k);
for ii = 1:K
    W = exprnd(1);
    X = mvnrnd(zeros(1,k), Sigma_(:,:,ii));
    vLaplace_rnd(ii,:) = sqrt(W)*X;
end

end
