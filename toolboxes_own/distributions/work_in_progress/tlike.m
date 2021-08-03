function [nLogL,lls]=tlike(x,mu,sigma_,df_t)
% Negative Log likelihood of the Standardized T distribution
%
% USAGE:
%   [LL,LLS]=tlike(X,MU,VAR_,NU)
%
% INPUTS:
%   X      - Standardized T random variables, either scalar or column vector
%   MU     - Mean of X, either scalar or size(x)
%   SIGMA_ - Variance Parameter of X, either scalar or size(x). (NOT THE
%            ACTUAL VARIANCE)
%   df_t   - Degree of freedom parameters, either scalar or size(x)
%
% OUTPUTS:
%   LL    - Log-likelihood evaluated at X
%   LLS   - Vector of log-likelihoods corresponding to X
%
% COMMENTS:
%   df_t>2
%
% REFERENCES:
%   [1] Cassella and Berger (1990) 'Statistical Inference'
%   [2] MATLAB: t Location-Scale Distribution
%
% See also STDTCDF, STDTINV, STDTRND, STDTPDF

% Copyright:
% Kevin Sheppard
% kevin.sheppard@economics.ox.ac.uk
% Revision: 1    Date: 9/1/2004
%
% Modified:
% Michael Stollenwerk
% michael.stollenwerk@live.com
% 29.04.2017


[T,K]=size(x);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input Checking
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if K~=1
    error('X must be a column vector');
end

if nargin==4
    if length(mu)~=1 && ~all(size(mu)==[T K])
        error('mu must be either a scalar or the same size as X');
    end
    if any(sigma_<=0)
        error('sigma2 must contain only positive elements')
    end
    if length(sigma_)==1
        sigma_=sigma_*ones(T,K);
    elseif size(sigma_,1)~=T || size(sigma_,2)~=1
        error('sigma2 must be a scalar or a vector with the same dimensions as X');
    end
    if length(df_t) == 1
        if df_t<=2
            error('NU must be greater than 2');
        end
    elseif size(df_t,1)~=T || size(df_t,2)~=1 || any(df_t <= 2)
        error('All values in NU must be greater than 2, and NU must be T by 1');
    end
    x=x-mu;
else
    error('Only 4 inputs supported');
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input Checking
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Compute the log likelihood
lls = gammaln(0.5.*(df_t+1)) - gammaln(df_t./2) - 1/2.*log(pi.*(df_t))...
    - 0.5.*(log(sigma_)) - ((df_t+1)./2).*(log(1 + (x.^2)./(sigma_.*(df_t))));
nLogL = -sum(lls);
