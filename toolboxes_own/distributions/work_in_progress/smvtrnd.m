function [ smvt_rnd ] = smvtrnd( df_t, Sigma_)
%SMVSTRND Generate STANDARDIZED t random array. (rows correspond to draws)
%
% USAGE:
%  [ smvst_rnd ] = smvstrnd( df_t, skew_, Sigma_)
%
% INPUTS:
%   DF_T         - N by 1 array of degree of freedom parameters.
%   SKEW_        - N by p array. Each row regulates the skewnness of the 
%                  multivariate skew t distribution. 
%   SIGMA_        - p by p by N array. Each p by p matrix corresponds to one 
%                  parameter matrix which regulates the dispersion matrix.
%
% OUTPUTS:
%   SMVT_RND    - N by p array. Each row corresponds to one random
%                  draw according to the respective ALPHA_, SIGMA_ and DF.
%
% See also MVNRND (Statistics and Machine Learning Toolbox - R2017b)
%
% REFERENCES:
%      [1] Stochastic Representation in our notes on the skew-t F project.

% Michael Stollenwerk
% michael.stollenwerk@live.com
% 21.07.2018

% DEPENDENCIES:
%   MVNRND (Statistics and Machine Learning Toolbox - R2017b)

%% Error checking
narginchk(2,2);

if numel(size(Sigma_)) ~=2 && numel(size(Sigma_)) ~=3
    error('Sigma dimensions are not valid')
end

if size(Sigma_,1) ~= size(Sigma_,2)
    error('Sigmas are not quadratic matrices')
end

if Sigma_ ~= permute(Sigma_, [2,1,3])
    error('Sigma Matrices are not symmetric')
end

[~,K,N] = size(Sigma_);

if any(df_t == inf) || any(df_t < K)
    error('All df_t must be K < df_t < inf')
end
%% Random Draws
mvn_rnd = mvnrnd(zeros(K,1), eye(K), N);
chi2_rnd = chi2rnd(df_t);
%% Stochastic Representation
smvt_rnd = NaN(N,K);

for ii = 1:N
    smvt_rnd(ii,:) = ...
        sqrt((df_t(ii) - 2)/df_t(ii)) / sqrt(chi2_rnd(ii)/df_t(ii)) * ...
            sqrtm(Sigma_(:,:,ii)) * mvn_rnd(ii,:)';
end  
