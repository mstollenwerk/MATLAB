function [ smvst_rnd ] = smvstrnd( df_t, skew_, Sigma_)
%SMVSTRND Generate STANDARDIZED skew-t random array. (rows correspond to
%draws)
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
%   SMVST_RND    - N by p array. Each row corresponds to one random
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
narginchk(3,3);
if size(skew_,2) == 1 && (size(skew_,1) > 1)
    skew_ = skew_';
end

if numel(size(Sigma_)) ~=2 && numel(size(Sigma_)) ~=3
    error('Sigma dimensions are not valid')
end

if size(Sigma_,1) ~= size(Sigma_,2)
    error('Sigmas are not quadratic matrices')
end

if Sigma_ ~= permute(Sigma_, [2,1,3])
    error('Sigma Matrices are not symmetric')
end

[~,p,N] = size(Sigma_);

if any(df_t == inf) || any(df_t < p)
    error('All df_t must be p < df_t < inf')
end
%% Random Draws
mvn_rnd = mvnrnd(zeros(p,1), eye(p), N);
sn_rnd = randn(N,1);
chi2_rnd = chi2rnd(df_t);
%% Stochastic Representation
smvst_rnd = NaN(N,p);

for ii = 1:N
    smvst_rnd(ii,:) = ...
        sqrt((df_t(ii) - 2)/df_t(ii)) / sqrt(chi2_rnd(ii)/df_t(ii)) * ...
        (...
            sqrtm( inv( Sigma_(:,:,ii) ) + skew_(ii,:)'*skew_(ii,:) ) \...
            mvn_rnd(ii,:)' +...
            Sigma_(:,:,ii)*skew_(ii,:)' /...
            sqrt(1 + skew_(ii,:)*Sigma_(:,:,ii)*skew_(ii,:)') *...
            abs(sn_rnd(ii))...
        );
end  
