function [ nLogL, logLcontr ] = smvstlike( data, df_t, Sigma_, alpha_ )
%SMVSTLIKE Negative log-likelihood for the STANDARDIZED multivariate skew-t.
%
% USAGE:
%   [ nLogL, logLcontr ] = smvstlike( data, df_t, Sigma_, alpha_ )
%
% INPUTS:
%   DATA         - N by p array of vector realizations.
%   DF           - N by 1 array of degree of freedom parameters.
%   ALPHA_       - N by p array. Each row regulates the skewnness of the 
%                  multivariate skew t distribution. 
%   SIGMA_       - p by p by N array. Each p by p matrix corresponds to one 
%                  parameter matrix which regulates the covariance matrix.
%
% OUTPUTS:
%   NLOGL        - Negative log-likelihood value
%   LOGLCONTR    - Log-likelihood contributions
%
% COMMENTS:
%
% REFERENCES:
%      [1] Bodnara, Conrad, Parolya and Stollenwerk (2018)
%      [2] my_notes
%
% DEPENDENCIES:
%  tcdf [Statistics and Machine Learning Toolbox (2017b)]
%
%  See also mvstlike, tcdf [Statistics and Machine Learning Toolbox (2017b)]
%
% Michael Stollenwerk
% michael.stollenwerk@live.com
% 19.07.2018
% 21.07.2018: - switched 2d-array dimensions to match MATLAB syntax
%             - adjusted for change of mvtlike.
% 17.08.2018: - changed documentation

% warning('smvstlike has not been properly tested for bugs yet')

%% Error checking
narginchk(4,4);

if numel(size(Sigma_)) ~=2 && numel(size(Sigma_)) ~=3
    error('Sigma dimensions are not valid')
end

if size(Sigma_,1) ~= size(Sigma_,2)
    error('Sigmas are not quadratic matrices')
end

if Sigma_ ~= permute(Sigma_, [2,1,3])
    error('Sigma Matrices are not symmetric')
end

if (size(Sigma_,1) ~= size(data,2)) || (size(Sigma_,2) ~= size(data,2)) || (size(data,1) ~= size(alpha_,1)) || (size(data,2) ~= size(alpha_,2))
   error('Size of Sigma, data and alpha do not conform')
end

[~,p,N] = size(Sigma_);
%% Log-likelihood computation

logLcontr=NaN(N,1);

for ii = 1:N
    
    nlogLcontr_mvt = mvtlike( ...
        [],Sigma_(:,:,ii)*(df_t(ii) - 2)/df_t(ii),df_t(ii),data(ii,:)' ...
    );

    x = alpha_(ii,:) * data(ii,:)' / ...
        sqrt(df_t(ii) - 2 + data(ii,:)/Sigma_(:,:,ii)*data(ii,:)') *... 
        sqrt(df_t(ii) + p);
    
    logLcontr(ii) = log(2) - nlogLcontr_mvt + log(tcdf(x, df_t(ii) + p));
    
end

nLogL = -sum(logLcontr);
end