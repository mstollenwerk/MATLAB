function [ ncWish_rnd ] = ncWishrnd_( df, Sigma_, Omega_)
%MNCFLIKE Negative log-likelihood for the matrix-variate non-central F distr.
%
% USAGE:
%  [ mncf_rnd ] = ncWishrnd( nu1, Sigma_, Omega_)
%
% INPUTS:
%   DF_1         - Vector of length N. First degrees of freedom parameters 
%                  of rescaled matrix variate non-central F distribution
%   SIGMA_       - p by p by N array of parameter matrices, 
%                  regulates the covariance. Matrix Parameter present in 
%                  both noncentral F and skew t distribution.
%   OMEGA_       - p by p by N array of non-centrality matrices of rescaled
%                  matrix variate non-central F distribution
%
% OUTPUTS:
%   MNCF_RND     - Random Matrix
%
% COMMENTS:
%   See also iwishrnd, mvnrnd [Statistics and Machine Learning Toolbox R2017b]
%
% REFERENCES:
%      [1] 
% DEPENDENCIES:
%   iwishrnd, mvnrnd [Statistics and Machine Learning Toolbox R2017b]

% Michael Stollenwerk
% michael.stollenwerk@live.com
% 21.07.2018

%% Input Checking
% add Sigma checks
if numel(size(Omega_)) ~= 2 && numel(size(Omega_)) ~= 3
    error('Omega_ dimensions are not valid')
end

if size(Omega_,1) ~= size(Omega_,2)
    error('Omega_ are not a quadratic matrices')
end

if Omega_ ~= permute(Omega_, [2,1,3])
    error('Omega Matrices are not symmetric')
end

[p,~,N] = size(Omega_);

for ii=1:N
    try
        chol(Sigma_(:,:,ii));
    catch
        error('Sigma Matrices are not positive definite')
    end
end

if any(mod(df,1) ~= 0)
%     warning('Rounded df_1')
    df = round(df);
end
%% Random Draws
ncWish_rnd = NaN(p,p,N);

for ii=1:N
    M_sq  = Sigma_(:,:,ii)*Omega_(:,:,ii);
    M = [sqrtm(M_sq), zeros(p,df(ii)-p)];

    X = reshape(mvnrnd(M(:), kron(Sigma_(:,:,ii),eye(df(ii)))), p, df(ii));
            
    ncWish_rnd(:,:,ii) = X*X';
end
end
