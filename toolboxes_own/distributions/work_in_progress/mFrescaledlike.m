function [ nLogL, logLcontr ] = mFrescaledlike( data_mat, df_1, df_2_rescaled, Sigma_rescaled)
%MNCFLIKE Negative log-likelihood for the matrix-variate F distr.
%
% USAGE:
%  [ nLogL, logLcontr ] = mFlike( data_mat, df_1, df_2, Sigma_)
%
% INPUTS:
%   data_mat      - p by p by N array of realized covariance matrix data
%   DF_1          - Vector of length N. First degrees of freedom parameters 
%                   of matrix variate F distribution.
%   DF_2          - Vector of length N. Second degrees of freedom parameter
%                   of matrix variate F distribution.
%   SIGMA_        - p by p by N array of parameter matrices, 
%                   regulates the covariance. 
%
% OUTPUTS:
%   NLOGL        - Negative log-likelihood value
%   LOGLCONTR    - Log-likelihood contributions
%
% See also MVGAMMALN MFRND FMEST
%
% COMMENTS:
%
% REFERENCES:
%   [1] Gupta and Nagar (2001) - Matrix Variate Distributions, p.156.
%   [2] Mulder and Pericchi (2018) - The Matrix F Prior for Estimating and
%           Testing Covariance Matrices.
%   [3] Opschoor, Janus, Lucas and Van Dijk (2018) - New HEAVY Models for 
%           Fat-Tailed Realized Covariances and Returns.
%
% DEPENDENCIES:
%   MVGAMMALN
%
% Michael Stollenwerk
% michael.stollenwerk@live.com
% 25.02.2020
%% Input Checking
if numel(size(data_mat)) ~=2 && numel(size(data_mat)) ~=3
    error('data_mat dimensions are not valid')
end
if numel(size(Sigma_rescaled)) ~=2 && numel(size(Sigma_rescaled)) ~=3
    error('Sigma_ dimensions are not valid')
end

if size(data_mat,1) ~= size(data_mat,2)
    error('data_mat are not a quadratic matrices')
end
if size(Sigma_rescaled,1) ~= size(Sigma_rescaled,2)
    error('Sigma_ are not a quadratic matrices')
end

if data_mat ~= permute(data_mat, [2,1,3])
    error('data_mat Matrices are not symmetric')
end
if Sigma_rescaled ~= permute(Sigma_rescaled, [2,1,3])
    error('Sigma_ Matrices are not symmetric')
end

if size(Sigma_rescaled) ~= size(data_mat)
    error('Size of Sigma and data_mat does not match')
end

[p,~,N] = size(data_mat);
%% Log-likelihood computation
logLcontr=NaN(N,1);

for ii = 1:N

    A = data_mat(:,:,ii);
    Sig = Sigma_rescaled(:,:,ii);
    nu1 = df_1(ii);
    nu2 = df_2_rescaled(ii);
    invSig_c = nu1 / (nu2 - p - 1) * inv(Sig);
    
    term1 = -mvbetaln(nu1/2, nu2/2, p);
    term2 = nu1 / 2 * log(det(invSig_c));
    log_normalizing_constant = term1 + term2;
      
    term3 = (nu1 - p - 1)/2*log(det(A));
    term4 = -(nu1 + nu2)/2*log(det(eye(p) + invSig_c * A));
    log_kernel = term3 + term4;

    logLcontr(ii) = log_normalizing_constant + log_kernel;    
    
end

nLogL = -sum(logLcontr);
