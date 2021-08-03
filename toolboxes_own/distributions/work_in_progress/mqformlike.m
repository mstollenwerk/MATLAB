function [ nLogL, logLcontr ] = mqformlike( data_mat, Sigma_, Psi_, A, mhg_precision)
%MQFORMLIKE Negative log-likelihood for the distribution of matrix quadratic 
%forms of the normal distribution.
%
% USAGE:
%  [ nLogL, logLcontr ] = mqformlike( data_mat, Sigma_, Psi_, A, mhg_precision)
%
% INPUTS:
%   data_mat      - p by p by N array of realized covariance matrix data.
%   SIGMA_        - p by p by N array of parameter matrices.
%   PSI_          - df_ by df_ by N array of parameter matrices.
%   A             - df_ by df_ by N array of paramter matrices.
%   MHG_PRECISION - Number of Jack Functions used to approximate the matrix
%                   valued hypergeometric function. This input is inversly
%                   related to computing time.
%
% FROM HERE ON DOCUMENTATION HAS TO BE UPDATED!!!!!!!!  
%
% OUTPUTS:
%   NLOGL        - Negative log-likelihood value
%   LOGLCONTR    - Log-likelihood contributions
%
% See also MVGAMMALN LOGMHG.C LOGMHG
%
% COMMENTS:
%   Description of LOGMHG:
%       Syntax: 	[s,c]=mhg([M K lambda],alpha,a,b,x,y)
%       Description: 	Computes the truncated Hypergeometric function of one or two matrix arguments
%       pFq(alpha)( a1,..., ap; b1,..., bq; X, Y) as a series of Jack functions, truncated for partitions of size not exceeding M.
%       Arguments: 	
% 
%           alpha is positive real
%            (typically alpha=2 in settings involving REAL random matrices, and alpha=1 in settings involving COMPLEX random matrices);
%           a and b are real arrays;
%           x and y are real arrays containing the eigenvalues of the matrix arguments X and Y, respectively;
%           The argument y may be omitted.
%           The argument K may be omitted but if specified, the summation is over partitions kappa such that kappa[1]<=K.
%           The argument lambda may be omitted, but if specified, the summation is over partitions kappa such that kappa[i]<=lambda[i] for i=1,2,... 
% 
%       Output: 	
% 
%           s is the value of the truncated hypergeometric function
%           c is a vector of size M+1 with the marginal sums for partitions of size 0,1,...,M (thus s=sum(c)). 
% 
%       Comments: 	Larger values of M will yield more accurate results. There are no rules on selecting the optimal M. Start with (say) M=10 and experiment until you obtain the best value of M for your application. The values of the output vector c may give you some idea about the convergence.
%       Example 1: 	mhg(20,2,[0.1,0.2],[0.3,0.4],[0.5 0.6 0.7]) returns 3.1349, which is a good approximation to
%       2F2(2)(0.1,0.2; 0.3,0.4; diag(0.5,0.6,0.7)).
%       Example 2: 	mhg(20,1,[ ],[0.1],[0.2 0.3 0.4],[0.5 0.6 0.7]) returns 6.6265, which is a good approximation to
%       0F1(1)(0.1; diag(0.2,0.3,0.4),diag(0.5,0.6,0.7)).
%       Example 3: 	mhg([20 4],1,[0.8],[0.1],[0.2 0.3 0.4]) returns 13.4183, which is an approximation to
%       1F1(1)(0.8;0.1; diag(0.2,0.3,0.4)) summed over partitions whose size does not exceed 20 and all of whose parts do not exceed 4.
%       Example 4: 	mhg([10 4 4 3 1],1,[0.8],[0.1],[0.2 0.3 0.4]) returns 13.4182, which is an approximation to
%       1F1(1)(0.8;0.1; diag(0.2,0.3,0.4)) summed over partitions kappa whose size does not exceed 20 AND all of whose parts do not exceed 4 AND kappa[1]<=4, kappa[2]<=3, kappa[3]<=1. 
%
% REFERENCES:
%      [1] Bodnara, Conrad, Parolya and Stollenwerk (2018)
%      [2] Koev and Edelman (2006), The Efficient Evaluation of the 
%          Hypergeometric Function of a Matrix Argument. 
%          Code:
%          http://www.math.sjsu.edu/~koev/software/mhgref.html

% Michael Stollenwerk
% michael.stollenwerk@live.com
% 02.12.2019

% DEPENDENCIES:
%   MVGAMMALN LOGMHG.C LOGMHG

[k,~,N] = size(Sigma_);
[df_, ~, ~] = size(A);

%% Check logmhg Installation.
% Code by Koev and Edelman (2006)
try
    log(mhg(1,2,[0.1,0.2],[0.3,0.4],[0.5 0.6 0.7]));
catch
    disp(' ')
    disp('==============================================================================================')
    disp(' The function MHG (Hypergeometric Function of a Matrix Argument) has to be installed.')
    disp(' ')
    disp(' 1. Make sure you have a c++ compiler installed on your system. ')
    disp(' 2. Download and unzip http://www.math.sjsu.edu/~koev/software/mhg15.zip onto your Matlab path. ')
    disp(' 3. Run "mex mhg.c" at the MATLAB prompt. ')
    disp('==============================================================================================')
    error('logmhg problem. Siehe oben')
end
%% Input Checking
% Are matrix dimensions valid?
if numel(size(data_mat)) ~=2 && numel(size(data_mat)) ~=3
    error('data_mat dimensions are not valid')
end
if numel(size(Sigma_)) ~=2 && numel(size(Sigma_)) ~=3
    error('Sigma_ dimensions are not valid')
end
if numel(size(Psi_)) ~=2 && numel(size(Psi_)) ~=3
    error('Psi_ dimensions are not valid')
end
if numel(size(A)) ~=2 && numel(size(A)) ~=3
    error('A dimensions are not valid')
end

% Are matrices quadratic?
if size(data_mat,1) ~= size(data_mat,2)
    error('data_mat are not a quadratic matrices')
end
if size(Sigma_,1) ~= size(Sigma_,2)
    error('Sigma_ are not a quadratic matrices')
end
if size(Psi_,1) ~= size(Psi_,2)
    error('Psi_ are not a quadratic matrices')
end

% Are respective matrices symmetric?
if data_mat ~= permute(data_mat, [2,1,3])
    error('data_mat Matrices are not symmetric')
end
if Sigma_ ~= permute(Sigma_, [2,1,3])
    error('Sigma_ Matrices are not symmetric')
end
if Psi_ ~= permute(Psi_, [2,1,3])
    error('Psi_ Matrices are not symmetric')
end

% Do matrix dimensions match ?
if size(Sigma_) ~= size(data_mat)
    error('Size of Sigma_ and data_mat does not match')
end
if size(Psi_) ~= size(data_mat)
    error('Size of Psi_ and data_mat does not match')
end
if size(Psi_) ~= size(Sigma_)
    error('Size of Omega_ and Sigma_ does not match')
end

if isempty(mhg_precision)
    mhg_precision = 25;
end

% disp('MHG precision parameter (number of jack functions to be)')
% disp(['calculated is set to ', num2str(mhg_precision), '. This is strongly related to'])
% disp('computing time.')

%% Log-likelihood computation
logLcontr=NaN(N,1);

for ii = 1:N
    
    inv_sqrtA = inv(sqrtm(A(:,:,ii)));
    B = eye(df_) - inv_sqrtA*inv(Psi_(:,:,ii))*inv_sqrtA;
    
    logLcontr(ii) = -.5*df_*k*log(2) - mvgammaln(.5*df_, k) ...
        - .5*k*log(det(A(:,:,ii)*Psi_(:,:,ii))) - .5*df_*log(det(Sigma_(:,:,ii))) ...
        + .5*(df_ - k - 1)*log(det(data_mat(:,:,ii))) + trace(-.5*(Sigma_(:,:,ii)\data_mat(:,:,ii))) ...
        + log(mhg(mhg_precision, df_, [], [], eig(B), eig(.5*Sigma_(:,:,ii)\data_mat(:,:,ii))));    
end

nLogL = -sum(logLcontr);
