function [ weights ] = gmvp_noshort( covmat )
%GMVP_NOSHORT calculates the globl minimum variance portfolio weights without shorting.
%   Detailed explanation goes here
%
% DEPENDENCIES: 
%   Financial Toolbox
%
% Michael Stollenwerk
% michael.stollenwerk@live.com
% 25.04.2019

[N,~,T] = size(covmat);

weights = NaN(T,N);
for t=1:T
    P = Portfolio('AssetCovar', covmat(:,:,t), 'LowerBound', 0, 'Budget', 1,...
        'AssetMean', zeros(N,1)); %AssetMean does not matter but has to be set.
    weights(t,:) = estimateFrontierLimits(P,'Min'); 
end

end
