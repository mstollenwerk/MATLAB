function [ nLogL, logLcontr, varargout ] = ...
    matvgWisht3like( Sigma_, Theta_, df_t, X, varargin )
%MATVGWISHTLIKE Negative log-likelihood, score and fisher information 
%   matrix of the matrix-variate generalized Wishart distribution based on 
%   multivariate t-distribution.
%
% USAGE:
%   [ nLogL, score, param ] = matvgWishtlike( Sigma_, Theta_, df_n, df_n, X )
%   [ nLogL, score, param ] = matvgWishtlike( [], [], [], [], X, all_param )
%
% INPUTS:
%   X           - Array (p by p by N). Symmetric p.d. data matrices.
%   SIGMA_      - Array (p by p). Symmetric p.d. parameter matrix.
%   THETA_      - Array (p by p). Symmetric p.d. parameter matrix.
%   DF_n        - Double. First degrees of freedom parameter. Corresponds
%                 to "number of sums".
%   DF_2        - Double. Second degrees of freedom parameter. Corresponds
%                 to degree of freedom parameter of underlying t-distribution.
%   ALL_PARAM   - 
%
% OUTPUTS:
%   NLOGL       - Double. Negative log-likelihood value.
%   LOGLCONTR   - Array (N by 1). Log-likelihood contriubtions.
%   SCORE       - Struct. Fields as parameter names. Contains derivatives
%                 of log-likelihood w.r.t. parameters.
%   FISHERINFO  - 
%   PARAM       -
%   
%
% See also 
% 
% COMMENTS:
%
%
% REFERENCES:
%   [1] Díaz-García (2013) - Distribution theory of quadratic forms for
%       matrix multivariate elliptical distributions.
%
% DEPENDENCIES:
%
%
% Michael Stollenwerk
% michael.stollenwerk@live.com
% 14.02.2021


[p,~,N] = size(X);
p_ = p*(p+1)/2;
narginchk(4,5); %%%%%%%
nargoutchk(0,5);
%% Param
if nargin == 5 %%%%%%%
    if ~(isempty(Sigma_) && isempty(Theta_) && isempty(df_t))
        error('Cannot input all_param and any of the parameters individually!')
    end
    all_param = varargin{1};
    Sigma_ = ivechchol(all_param(1 : p_));
    
    df_n_ = length(all_param) - p_ - 1;
    df_n = (2*df_n_ + 1/4)^(1/2) - 1/2;
    
    Theta_ = ivechchol(all_param(p_ + 1 : p_ + df_n_));
    df_t = all_param(p_ + df_n_ + 1);
end

% Optional Parameter Output
if nargout == 5
    param.Sigma_ = Sigma_;
    param.Theta_ = Theta_;
    param.chol_Sigma_ = vechchol(Sigma_);
    param.chol_Theta_ = vechchol(Theta_);
    param.df_t = df_t;
    param.all = [param.chol_Sigma_; param.chol_Theta_; df_t];
    
    varargout{3} = param;
end
%% Log-likelihood computation
logLcontr = NaN(N,1);

term1 = gammaln( (df_t + df_n*p)/2 ) ...
        - gammaln(df_t/2) ...
        - mvgammaln(df_n/2, p);
term2 = df_n*p/2*log(df_t);
term3 = -df_n/2*log(det(Sigma_));
term4 = -p/2*log(det(Theta_));

log_norm_const = term1 + term2 + term3 + term4;

for ii = 1:N
    A = X(:,:,ii);
    
    term5 = (df_n-p-1)/2*log(det(A));
        
    % Only the nonzero eigenvalues are taken (if not all -> spiked ncF).
        x_mhg_1 = eig(Theta_).^(-1);
        x_mhg_2 = eig((Sigma_\A)./df_t);
        
        try
            y_mhg = mhg(...
                9, 2, ...
                (df_n*p + df_t)/2, [], x_mhg_1, x_mhg_2 ...
            );
            term6 = log(y_mhg);
        catch ME
            disp(ME)
            disp('Error in MHG function evaluation. Is it installed?')
            disp(' ')
            disp('==============================================================================================')
            disp(' The function MHG (Hypergeometric Function of a Matrix Argument) has to be installed.')
            disp(' ')
            disp(' 1. Make sure you have a c++ compiler installed on your system. ')
            disp(' 2. Download and unzip http://www.math.sjsu.edu/~koev/software/mhg15.zip onto your Matlab path. ')
            disp(' 3. Run "mex mhg.c" at the MATLAB prompt. ')
            disp('==============================================================================================')
            error('mhg problem. Siehe oben')
        end
    
    log_kernel = term5 + term6;
		
    logLcontr(ii) = log_norm_const + log_kernel;
end
nLogL = -sum(logLcontr);
%% Score computation
if nargout >= 3
%     invSig = inv(Sigma_);
%     
%     scoreSig = NaN(N,p_);
%     score_df_1 = NaN(N,1);
%     score_df_2 = NaN(N,1);
%     for ii = 1:N
%         
%         % General matrix derivative (ignoring symmetry of Sigma_): [ enter to matrixcalculus.org: (df2 + p - 1)/2 * log(det(Sigma))-(df1 + df2 + p - 1)/2*log(det (Sigma + X )) ]
% 
%         % Apply vec operator.
%         S = S(:);
%         
%         % Accounting for symmetry of Sigma_:
%         D = Dmatrix(p);
%         scoreSig(ii,:) = D'*S;
% 
%         % Score nu
%         score_df_1(ii) = NaN;
%         score_df_2(ii) = NaN;
% 
%     end
%     
%     score.Sigma_ = scoreSig;
%     score.df_1 = score_df_1;
%     score.df_2 = score_df_2;
%     
%     varargout{1} = score;
varargout{1} = NaN;

end
%% Fisher Info
if nargout >= 4
%     
%     fisherinfo.Sigma_ = NaN(p_,p_,N);
%     fisherinfo.df_1 = NaN(N,1);
%     fisherinfo.df_2 = NaN(N,1);
%     
%     varargout{2} = fisherinfo;
varargout{2} = NaN;
end

end
