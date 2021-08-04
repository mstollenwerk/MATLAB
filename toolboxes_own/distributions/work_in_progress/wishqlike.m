function [ nLogL, logLcontr ] = wishqlike( sigma, data )
%WISHLIKE Negative quasi log-likelihood for the Wishart distribution.
%
% USAGE:
%   NLOGL = wishlike(sigma,df,data)
%
% INPUTS:
%   SIGMA        - K by K by N array of scale matrices of the Wishart
%                  distribution.
%   DATA         - K by K by N array of matrix data realizations.
%
% OUTPUTS:
%   NLOGL        - Negative log-likelihood value
%   LOGLCONTR    - Log-likelihood contributions
%
% COMMENTS:
%
%  See also 
%
% REFERENCES:
%      [1] Sutradhar and Ali (1989) - A Generalization of the Wishart
%      Distribution for the Elliptical Model and its Moments for the
%      Multivariate t Model.

% Michael Stollenwerk
% michael.stollenwerk@live.com
% 21.06.2019
%% Error checking
narginchk(2,2);

[k,k2,n] = size(sigma);
if k~=k2
    error('Covariance matrix parameters of Wishart distribution must be square, symmetric, and positive semi-definite.')
end

[k,k2,~] = size(data);
if k~=k2
    error('Covariance matrix data of Wishart distribution must be square, symmetric, and positive semi-definite.')
end

if size(sigma)~=size(data)
    error('Covariance matrix data and parameters of Wishart distribution must have the same dimension size')
end
%% Log-likelihood computation

logLcontr=NaN(n,1);

if isnumeric(sigma) && isnumeric(data) % Double input
    for i = 1:n 
        logLcontr(i) = -.5*log(det(sigma(:,:,i))) ...
                       -.5*trace(sigma(:,:,i)\data(:,:,i));
    end
elseif iscell(sigma) && iscell(data) % Cell input
    for i = 1:n
        logLcontr(i) = -.5*log(det(sigma{i})) ...
                       -.5*trace(sigma{i}\data{i});
    end
end

nLogL=-sum(logLcontr);
end

