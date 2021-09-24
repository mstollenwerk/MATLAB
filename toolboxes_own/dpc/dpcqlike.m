function [ nLogL, logLcontr ] = dpcqlike( d, diagLRL )
%DPCQLIKE calculates negative quasi log-likelihood of DPC model.
%
% USAGE:
%
% INPUTS:
%   d            - A N by T matrix of conditional components
%   R            - A N by T matrix of diag(L_t'*R_t*L_t), where L are the
%                  conditional loadings and R the covariance matrix
%
% OUTPUTS:
%   NLOGL        - Negative of log-likelihood
%   LOGLCONTR    - log-likelihood contributions
%
% COMMENTS:
%
%  See also DPCLIKE, SDPCREC
%
% REFERENCES:
%      [1] Aielli, G. P. and M. Caporin (2015) - Dynamic Principal
%      Components: a New Class of Multivariate GARCH Models

% Michael Stollenwerk
% michael.stollenwerk@live.com
% 02.09.2016
%% Input Checks
if size(d)~=size(diagLRL) % Both inputs must be same size
    error('d and diagLRL have to have the same size')
end
if max(size(d))~=max(size(diagLRL)) % All inputs must have same length T.
    error('All inputs must have same length')
end
%%
logLcontr=-.5*sum(log(d)+diagLRL./d);
logLcontr=logLcontr';
nLogL=-sum(logLcontr);
end
