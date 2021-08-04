function [ Sigma_, Psi_, A, nLogL ] = mqformest( data_mat, varargin )
%MQFORMEST Estimates the paramters of the distribution of matrix quadratic 
%forms of the normal distribution.
%
% USAGE:
%  [ Sigma_, Psi_, A, nLogL ] = mqformest( data_mat, varargin )
%
% INPUTS:
%   DATA         - p by p by N array of covariance matrix data.
%   VARARGIN     - 1. x0            - Optimization Starting Values
%                  2. mhg_precision - Number of Jack Functions used to 
%                                     approximate the matrix valued 
%                                     hypergeometric function. This input 
%                                     is inversly related to computing 
%                                     time.
%
%
% OUTPUTS:
%   SIGMA_       - p by p by N array of parameter matrices.
%   PSI_         - df_ by df_ by N array of paramter matrices.
%   A            - df_ by df_ by N array of paramter matrices.
%
% See also MQFORMLIKE MQFORMRND
%
% REFERENCES:

% Michael Stollenwerk
% michael.stollenwerk@live.com
% 19.12.2019

% DEPENDENCIES:
%   MQFORMLIKE

[k,~,~] = size(data_mat);
if isempty(varargin)
    mhg_precision = 25;
else
    mhg_precision = varargin{1};
end
%% 
Psi_x0 = eye(k);
Ax0 = eye(k);
Sigma_x0 = mean(data_mat,3)/k;
x0 = [vech(chol(Sigma_x0, 'lower'))', vech(chol(Psi_x0, 'lower'))', vech(chol(Ax0, 'lower'))'];
%%
eparam = fminunc( ...
    @(param) like_wrapper(param, data_mat, mhg_precision), x0, ...
    optimoptions('fminunc', 'Display', 'iter-detailed') ...
);
[nLogL, ~, Sigma_, Psi_, A ] = like_wrapper(eparam, data_mat, mhg_precision);
end
%%
function [nLogL, logLcontr, Sigma_, Psi_, A ] = like_wrapper(param, data_mat, mhg_precision) 

[k,~,T] = size(data_mat);
k_ = k*(k+1)/2;

Sigma_ = ivech(param(1:k_), 'lower');
Sigma_ = Sigma_ * Sigma_';
Psi_ = ivech(param(1+k_:k_+k_), 'lower');
Psi_ = Psi_ * Psi_';
A = ivech(param(1+k_+k_:k_+k_+k_), 'lower');
A = A*A';

[ nLogL, logLcontr ] = mqformlike( ... 
    data_mat, ... 
    repmat(Sigma_,1,1,T), ...
    repmat(Psi_,1,1,T), ...
    repmat(A,1,1,T), ...
    mhg_precision ...
);

end
