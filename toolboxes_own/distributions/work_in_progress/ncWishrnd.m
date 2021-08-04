function [ ncWish_rnd ] = ncWishrnd( Sigma_, Omega_, df_)
%NCWISHRND Generate non-central Wishart random matrices
%
% USAGE:
%  [ ncWish_rnd ] = ncWishrnd( df_, Sigma_, Omega_)
%
% INPUTS:
%   SIGMA_       - p by p by N array of parameter matrices, 
%                  regulates the covariance. 
%   OMEGA_       - p by p by N array of parameter matrices,
%                  regulates the non-centrality matrix.
%   DF_          - Vector of length N. Degree of freedom parameters 
%                  of non-central Wishart distribution.
%
% OUTPUTS:
%   NCWISH_RND   - Random Matrix
%
% COMMENTS:
%   OMEGA_ is NOT equal to the THETA matrix in Gupta, Nagar (1999), but 
%   THETA_ = inv(SIGMA_)*OMEGA_.
%
% REFERENCES:
%   [1] Gupta, Nagar (1999) - Matrix Variate Distributions, p.114.
%
% DEPENDENCIES:
%
%   See also 

% Michael Stollenwerk
% michael.stollenwerk@live.com
% 21.02.2019

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
% if any(df_ == inf) || any(df_ < p)
%     error('df_ must be p? < df_2 < inf')
% end

if any(mod(df_,1) ~= 0)
	warning('Rounded df_')
    df_ = round(df_);
end
%% Random Draws
ncWish_rnd = NaN(p,p,N);

% Copied from rwishart (lines 61-72) in R package matrixsampling.
for ii=1:N
    Sigma_root = sqrtm(Sigma_(:,:,ii));
    Omega_root = sqrtm(Omega_(:,:,ii));

    Z = randn(p);
    Y = randn(p, df_(ii) - p);
    
    ncWish_ = (Omega_root + Sigma_root*Z)*(Omega_root + Z'*Sigma_root) + ...
        Sigma_root*(Y*Y')*Sigma_root;
    
    if max(max(abs(ncWish_-ncWish_'))) < 1e-10
        ncWish_rnd(:,:,ii) = .5*(ncWish_+ncWish_'); % Force exact Symmetry.
    else
        ncWish_rnd(:,:,ii) = ncWish_;        
    end
end
end
