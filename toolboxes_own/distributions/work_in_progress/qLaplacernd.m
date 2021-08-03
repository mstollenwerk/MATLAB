function [ qLaplace_rnd ] = qLaplacernd( df, Sigma_ )
%QLAPLACERND
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
% 22.02.2020
%
% DEPENDENCIES:
%
%% 
[p, ~, k] = size(Sigma_);
[k_,~] = size(df);

if k~=k_
    error('Input array dimensions do not conform.')
end

qLaplace_rnd = NaN(p,p,k);
for ii = 1:k
    W = exprnd(1);
    Sigma_vecX = kron(Sigma_(:,:,ii), eye(df(ii)));
    vecX = mvnrnd(zeros(p*df(ii),1), Sigma_vecX);
    X = reshape(sqrt(W)*vecX,df(ii),p)';
    qLaplace_rnd(:,:,ii) = X*X';
end

end
