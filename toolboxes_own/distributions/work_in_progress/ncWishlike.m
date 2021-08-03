function [ nLogL, logLcontr ] = ncWishlike( data_mat, df, Sigma_, Omega_, mhg_precision)
%NCWISHRND Generate non-central Wishart random matrices
%
% USAGE:
%  [ ncWish_like ] = ncWishlike( data_mat, df, Sigma_, Omega_)
%
% INPUTS:
%   data_mat      - p by p by N array of realized covariance matrix data
%   DF_           - Vector of length N. Degrees of freedom parameters 
%                   of non-central Wishart distribution.
%   SIGMA_        - p by p by N array of parameter matrices, 
%                   regulates the covariance. 
%   OMEGA_        - p by p by N array of non-centrality matrices.
%
% OUTPUTS:
%   NLOGL        - Negative log-likelihood value
%   LOGLCONTR    - Log-likelihood contributions
%
% COMMENTS:
%   OMEGA_ is NOT equal to the THETA matrix in Gupta, Nagar (1999), but 
%   THETA_ = inv(SIGMA_)*OMEGA_.
%
% REFERENCES:
%   [1] Gupta, Nagar (1999) - Matrix Variate Distributions, p. 114.
%
% DEPENDENCIES:
%
%   See also NCWISHRND NCWISHEST

% Michael Stollenwerk
% michael.stollenwerk@live.com
% 14.03.2019

%% Input Checking
if numel(size(Sigma_)) ~=2 && numel(size(Sigma_)) ~=3
    error('Sigma_ dimensions are not valid')
end
if numel(size(Omega_)) ~=2 && numel(size(Omega_)) ~=3
    error('Omega_ dimensions are not valid')
end

if size(Sigma_,1) ~= size(Sigma_,2)
    error('Sigma_ are not a quadratic matrices')
end
if size(Omega_,1) ~= size(Omega_,2)
    error('Omega_ are not a quadratic matrices')
end

if Sigma_ ~= permute(Sigma_, [2,1,3])
    error('Sigma_ Matrices are not symmetric')
end
if Omega_ ~= permute(Omega_, [2,1,3])
    error('Omega_ Matrices are not symmetric')
end

if size(Omega_) ~= size(Sigma_)
    error('Size of Omega_ and Sigma_ does not match')
end

if Sigma_ ~= permute(Sigma_, [2,1,3])
    error('Omega Matrices are not symmetric')
end
if Omega_ ~= permute(Omega_, [2,1,3])
    error('Omega Matrices are not symmetric')
end

[p,~,N] = size(Omega_);

% for ii=1:N
%     try
%         chol(Sigma_(:,:,ii));
%     catch
%         error('Sigma Matrices are not positive definite')
%     end
% end
% 
% for ii=1:N
%     try
%         chol(Omega_(:,:,ii));
%     catch
%         error('Omega Matrices are not positive definite')
%     end
% end
% 
% if any(df == inf) || any(df < p)
%     error('df_ must be p? < df_2 < inf')
% end

if isempty(mhg_precision)
    mhg_precision = 100;
end

[p,~,N] = size(Omega_);
%% Log-likelihood computation
logLcontr=NaN(N,1);

for ii = 1:N
    Theta_ = Sigma_(:,:,ii)\Omega_(:,:,ii);
    X = Sigma_(:,:,ii)\data_mat(:,:,ii);
    logLcontr(ii) = ...
        - .5*p*df(ii)*log(2) ...
        - .5*df(ii)*log(det(Sigma_(:,:,ii))) ...
        - .5*trace(Theta_) ...
        - mvgammaln(df(ii)/2, p) ...
        + .5*(df(ii) - p - 1)*log(det(data_mat(:,:,ii))) ...
        - .5*trace(X) ...
        + log(mhg(mhg_precision, 2, [], .5*df(ii), eig(.25*Theta_*X)));
end

nLogL = -sum(logLcontr);
end
