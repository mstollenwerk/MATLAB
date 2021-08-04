function [ eparam ] = maLest_mm( data_mat )
%MALRND Generate multivariate-asymmetric Laplace random sample
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
eparam.mu_ = mean(data_mat);
eparam.Sigma_ = cov(data_mat) - eparam.mu_'*eparam.mu_;
end
