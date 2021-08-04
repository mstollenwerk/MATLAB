function [ maL_rnd ] = vaLaplacernd( mu_, Sigma_ )
%MALRND Generate vector-valued asymmetric Laplace random sample
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
if size(mu_,1) ~= size(Sigma_,3)
    error("mu and Sigma input size does not go together.")
end
if size(mu_,2) ~= size(Sigma_,1)
    error("mu and Sigma input size does not go together.")
end
[N, k] = size(mu_);

maL_rnd = NaN(N,k);
for ii = 1:N
    W = exprnd(1);
    X = mvnrnd(zeros(1,k), Sigma_(:,:,ii));
    maL_rnd(ii,:) = mu_(ii,:)*W + sqrt(W)*X;
end

end
