function [ nLogL, logLcontr, param, varargout ] = matvgbeta2like( Omega_, Psi_, a, b, X, varargin )
%
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
%      [1] Based on Gupta and Nagar p.167
%
% Michael Stollenwerk
% michael.stollenwerk@live.com
% 30.07.2020
%
% DEPENDENCIES:
%
%%
[p,~,T] = size(X);
p_ = p*(p+1)/2;
%% Parameters
if nargin == 6
    all_param = varargin{1};
    Omega_ = ivechchol(all_param(1:p_));
    Psi_ = ivechchol(all_param(p_+1:2*p_));
    a = all_param(2*p_ + 1);
    b = all_param(2*p_ + 2);
end

param.Omega_ = Omega_;
param.Psi_ = Psi_;
param.a = a;
param.b = b;
if all(Psi_==0,"all")
    param.all = [vech(chol(Omega_,'lower'))' zeros(1,p_) a b];
else
    param.all = [vech(chol(Omega_,'lower'))' vech(chol(Psi_,'lower'))' a b];
end

%% Log-Likelihood 
% normalizing constant
log_norm_constant = ...
    - mvbetaln(a,b,p) ...
    + b*log(det( Omega_ + Psi_ ));

% kernel
log_kernel = NaN(T,1);
for ii = 1:T
    
    Y = X(:,:,ii);
      
    log_kernel(ii) = ...
        ( a - .5*(p + 1) ) * log(det( Y - Psi_ )) ...
        - (a+b)*log(det( Omega_ + Y ));
    
end
logLcontr = log_norm_constant + log_kernel;
nLogL = -sum(logLcontr);
%% Score
if nargout == 4
    score.Omega_ = NaN(p);
    score.Psi_ = NaN(p);
    score.a = NaN;
    score.b = NaN;
    varargout = score;
end
end
