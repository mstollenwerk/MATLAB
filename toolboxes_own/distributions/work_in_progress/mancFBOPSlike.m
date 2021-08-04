function [ nLogL, score, param_trans ] = mancFBOPSlike( param, X, mhg_precision)
%MANCFBOPSLIKE Negative log-likelihood for the rescaled matrix-variate 
%   non-central F distribution as defined by 
%   Bodnar, Okhrin, Parolya and Stollenwerk (2020). 
%
% USAGE:
% 
%
% INPUTS:
%   PARAM
%     (1 : p(p+1)/2  ) - vech(chol(Sigma)): Sigma (p by p) regulates the covariance matrix.
%     (_ : p(p+1)    ) - vech(chol(Omega)): Omega (p by p) non-centrality matrix.
%     ( p(p+1)+1 )     - nu1: First degrees of freedom parameter.
%     ( p(p+1)+2 )     - nu2: Second degrees of freedom parameter.
%   X                  - p by p array, matrix realization.
%   MHG_PRECISION - Number of Jack Functions used to approximate the matrix
%                   valued hypergeometric function. This input is inversly
%                   related to computing time.
%
% OUTPUTS:
%   NLOGL        - Negative log-likelihood value
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
%      [1] Bodnar, Okhrin, Parolya and Stollenwerk (2020)
%      [2] Koev and Edelman (2006), The Efficient Evaluation of the 
%          Hypergeometric Function of a Matrix Argument. 
%          Code:
%          http://www.math.sjsu.edu/~koev/software/mhgref.html
%
% DEPENDENCIES:
%   MVGAMMALN LOGMHG.C LOGMHG
%
% See also MVGAMMALN LOGMHG.C LOGMHG
%
% Michael Stollenwerk
% michael.stollenwerk@live.com

p = size(X,1);
p_ = p*(p+1)/2;
%% Parameter Transformation
Sigma_ = ivech(param(1:p_), 'lower')* ...
         ivech(param(1:p_), 'lower')';
Omega_ = ivech(param(p_+1:2*p_), 'lower')* ...
         ivech(param(p_+1:2*p_), 'lower')';
nu1 = param(2*p_+1);
nu2 = param(2*p_+2);

param_trans.Sigma_ = Sigma_;
param_trans.Omega_ = Omega_;
param_trans.nu1 = nu1;
param_trans.nu2 = nu2;
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


if isempty(mhg_precision)
    mhg_precision = 25;
end

% disp('MHG precision parameter (number of jack functions to be)')
% disp(['calculated is set to ', num2str(mhg_precision), '. This is strongly related to'])
% disp('computing time.')

%% Log-likelihood computation
c_nu = nu1/(nu2 - p - 1);
Q = (Sigma_/X)*Sigma_ + c_nu*Sigma_;

logL = - mvbetaln(nu1/2, nu2/2, p) + p*nu1/2*log(c_nu) ...
    - (nu2+p+1)/2*log(det(X)) + (nu1+2*nu2)/2*log(det(Sigma_)) ...
    - (nu1+nu2)/2*log(det(Q)) - 1/2*trace(Sigma_\Omega_) ...
    + log(mhg(mhg_precision, 2, (nu1+nu2)/2, nu1/2, eig(1/2*c_nu*Omega_/Q)));

nLogL = -logL;
%% Score
% score = mancFBOPSscore(X, Sigma_, Omega_, nu1, nu2);
end