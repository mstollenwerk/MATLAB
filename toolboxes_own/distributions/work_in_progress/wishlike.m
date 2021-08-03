function [ nLogL, logLcontr ] = wishlike( sigma, df, data )
%WISHLIKE Negative log-likelihood for the Wishart distribution.
%
% USAGE:
%   NLOGL = wishlike(sigma,df,data)
%
% INPUTS:
%   SIGMA        - K by K by N array of scale matrices of the Wishart
%                  distribution
%   DF           - Degree of freedom of the Wishart distribution
%   DATA         - K by K by N array of covariance matrix realizations
%
% OUTPUTS:
%   NLOGL        - Negative log-likelihood value
%   LOGLCONTR    - Log-likelihood contributions
%
% COMMENTS:
%   This is not W(df, Sigma/df)!!!
%
%  See also 
%
% REFERENCES:
%      [1] Sutradhar and Ali (1989) - A Generalization of the Wishart
%      Distribution for the Elliptical Model and its Moments for the
%      Multivariate t Model.

% Michael Stollenwerk
% michael.stollenwerk@live.com
% 11.02.2017
%% Error checking
narginchk(3,3);

[k,k2,n] = size(sigma);
if k~=k2
    error('Covariance matrix parameters of t-Wishart distribution must be square, symmetric, and positive semi-definite.')
end

[k,k2,~] = size(data);
if k~=k2
    error('Covariance matrix data of t-Wishart distribution must be square, symmetric, and positive semi-definite.')
end

if size(sigma)~=size(data)
    error('Covariance matrix data and parameters of t-Wishart distribution must have the same dimension size')
end
%% Log-likelihood computation

logLcontr=NaN(n,1);

if isnumeric(sigma) && isnumeric(data) % Double input
    for i = 1:n 
        logLcontr(i) = -df/2*log(det(sigma(:,:,i))) - ...
                       1/2*trace(sigma(:,:,i)\data(:,:,i)) + ...
                       (df-k-1)/2*log(det(data(:,:,i)));
    end
elseif iscell(sigma) && iscell(data) % Cell input
    for i = 1:n
        logLcontr(i) = -df/2*log(det(sigma{i})) - ...
                       1/2*trace(sigma{i}\data{i}) + ...
                       (df-k-1)/2*log(det(data{i}));
    end
end

logLcontr = logLcontr - df*k/2*log(2) - mvgammaln(df/2,k);

nLogL=-sum(logLcontr);
end

