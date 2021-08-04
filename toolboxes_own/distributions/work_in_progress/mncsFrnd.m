function [ mncsF_rnd ] = mncsFrnd( df_1, df_2, Sigma_, Omega_)
%MNCFLIKE Negative log-likelihood for the matrix-variate non-central spiked F distr.
%
% USAGE:
%  [ mncf_rnd ] = mncFrnd( nu1, nu2, Sigma_, Omega_)
%
% INPUTS:
%   DF_1         - Vector of length N. First degrees of freedom parameters 
%                  of rescaled matrix variate non-central F distribution
%   DF_2         - Vector of length N. Second degrees of freedom parameter
%                  of rescaled matrix variate non-central F distribution
%   SIGMA_       - p by p by N array of parameter matrices, 
%                  regulates the covariance. Matrix Parameter present in 
%                  both noncentral F and skew t distribution.
%   OMEGA_       - p by p by N array of non-centrality matrices of rescaled
%                  matrix variate non-central F distribution
%
% OUTPUTS:
%   MNCSF_RND     - Random Matrix
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
% 05.03.2019

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

for ii=1:N
    try
        chol(Sigma_(:,:,ii));
    catch
        error('Sigma Matrices are not positive definite')
    end
end

if any(df_2 == inf) || any(df_2 < p)
    error('df_2 must be p? < df_2 < inf')
end

if any(mod(df_1,1) ~= 0)
%     warning('Rounded df_1')
    df_1 = round(df_1);
end
%% Random Draws
mncsF_rnd = NaN(p,p,N);

for ii=1:N
    if all(Omega_(:,:,ii)==0)
        Lambda_ = zeros(p,df_1(ii));
    else
        [sqrtOmega,res] = sqrtm(Omega_(:,:,ii));
        if max(max(imag(sqrtOmega))) < 1e-6
            sqrtOmega = real(sqrtOmega);
        else
            error('Omega square root is bad: Imaginary')
        end
        if isnan(res) || res > 1e-8
            error('Omega square root is bad')
        end
        Lambda_ = [sqrtOmega, zeros(p,df_1(ii)-p)];
    end
    sqrtm_Sigma = sqrtm(Sigma_(:,:,ii));
    inv_sqrtm_Sigma_Lambda = sqrtm_Sigma\Lambda_;
    X = reshape(mvnrnd(zeros(p*df_1(ii),1), eye(p*df_1(ii))), p, df_1(ii));
    V = (X + inv_sqrtm_Sigma_Lambda)*(X + inv_sqrtm_Sigma_Lambda)';
    sqrtm_V = sqrtm(V);
    invT = iwishrnd(eye(p), df_2(ii));
    mncsF_ = (df_2(ii) - p - 1)/df_1(ii) * sqrtm_Sigma*sqrtm_V*invT*sqrtm_V*sqrtm_Sigma;
    if max(max((abs((mncsF_ - mncsF_')./mncsF_)))) < 1e-08
        mncsF_ = (mncsF_ + mncsF_')/2;
%         warning([num2str(ii) ' is not exactly symmetric'])
    else
        error(['Not symmetric matrix produced in mncF_rnd MC rep ' num2str(ii)])
    end
    if any(isnan(mncsF_))
        error('NaNs produced')
    end
    mncsF_rnd(:,:,ii) = mncsF_;    
end
end
