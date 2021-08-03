function [ nLogL, logLcontr, varargout ] = ...
    matvFlike( Sigma_, df_1, df_2, X, varargin )
%MATVFLIKE Negative log-likelihood and score of the matrix-variate F distr.
%
% USAGE:
%   [ nLogL, score, param ] = matvFlike( Sigma_, df_1, df_2, X )
%   [ nLogL, score, param ] = matvFlike( [], [], [], X, all_param )
%
% INPUTS:
%   X       - Array (p by p). Symmetric p.d. data matrix.
%   SIGMA_  - Array (p by p). Symmetric p.d. parameter matrix, 
%             regulates covariance. 
%   DF_1    - Double. First degrees of freedom parameter. 
%   df_2    - Double. Second degrees of freedom parameter. 
%
% OUTPUTS:
%   NLOGL   - Double. Negative log-likelihood value.
%   SCORE   - Struct. Fields as parameter names. Contains derivatives
%             of log-likelihood w.r.t. parameters.
%
% See also MVBETALN MATVFRND IVECHCHOL
% 
% COMMENTS:
%
%   Derivatives were derived using matrixcalculus.org
%   References [3] and [4] use slightly different definitions of the
%   distributions.
%
% REFERENCES:
%   [1] Gupta and Nagar (2001) - Matrix Variate Distributions, p.167, p.156.
%   [2] Mulder and Pericchi (2018) - The Matrix F Prior for Estimating and
%           Testing Covariance Matrices.
%   [3] Opschoor, Janus, Lucas and Van Dijk (2018) - New HEAVY Models for 
%           Fat-Tailed Realized Covariances and Returns.
%   [4] Bodnar, Okhrin, Parolya and Stollenwerk (2020)
%
% DEPENDENCIES:
%   MVBETALN IVECHCHOL
%
% Michael Stollenwerk
% michael.stollenwerk@live.com
% 16.06.2020

[p,~,N] = size(X);
p_ = p*(p+1)/2;
narginchk(4,5);
nargoutchk(0,5);
%% Param
if nargin == 5
    if ~(isempty(Sigma_) && isempty(df_1) && isempty(df_2))
        error('Cannot input all_param and any of the parameters individually!')
    end
    all_param = varargin{1};
    Sigma_ = ivechchol(all_param(1 : p_));
    df_1 = all_param(p_ + 1);
    df_2 = all_param(p_ + 2);
end

% Optional Parameter Output
if nargout == 5
    param.Sigma_ = Sigma_;
    param.chol_Sigma_ = vechchol(Sigma_);
    param.df_1 = df_1;
    param.df_2 = df_2;
    param.all = [param.chol_Sigma_; df_1; df_2];
    
    varargout{3} = param;
end
%% Input Checking
cholSig = chol(Sigma_,'lower'); % Checking if Sigma_ is symmetric p.d.
%% Log-likelihood computation
logLcontr = NaN(N,1);
term1 = -mvbetaln(df_1/2, df_2/2, p);
term2 = df_2/2 * logdet(Sigma_);
log_normalizing_constant = term1 + term2;

for ii = 1:N
    term3 = (df_1 - p - 1)/2*logdet( X(:,:,ii) );
    term4 = -(df_1 + df_2)/2*logdet (Sigma_ + X(:,:,ii) );
    log_kernel = term3 + term4;

    logLcontr(ii) = log_normalizing_constant + log_kernel;
end
nLogL = -sum(logLcontr);
%% Score computation
if nargout >= 3
    invSig = inv(Sigma_);

    score.Sigma_ = NaN(N,p_);
    score.n = NaN(N,1);
    score.nu = NaN(N,1);
    
    for ii = 1:N
        
        % General matrix derivative (ignoring symmetry of Sigma_): [ enter to matrixcalculus.org: (df2 + p - 1)/2 * log(det(Sigma))-(df1 + df2 + p - 1)/2*log(det (Sigma + X )) ]
        S = df_2/2*invSig - (df_1 + df_2)/2*inv( Sigma_ + X(:,:,ii) );
%         S = - df_1/2*invSig ...
%             + (df_1+df_2)/2 ...
%               *invSig/(eye(p) + X(:,:,ii)*invSig)*X(:,:,ii)*invSig;

        % Accounting for symmetry of Sigma_:
        S = 2*S - diag(diag(S));
        
        % Numerical inaccuracies make S not exactly symmetric
        S = (S+S')./2;
                
        score.Sigma_(ii,:) = vech(S);

        % Score nu
        score.df_1(ii) = NaN;
        score.df_2(ii) = NaN;

    end
    
    varargout{1} = score;

end
%% Fisher Info
% enter to matrixcalculus.org: nu2/2*inv(Sigma) - (df1 + nu2)/2*inv( Sigma + X ), change the sign, then take expectation to arrive at
% take symmetry into account (matrixcalculus.org doesn't) [D' ____ D]
if nargout >= 4
    
    fisherinfo.Sigma_ = NaN;
    fisherinfo.df_1 = NaN;
    fisherinfo.df_2 = NaN;
    
    varargout{2} = fisherinfo;
end

end
