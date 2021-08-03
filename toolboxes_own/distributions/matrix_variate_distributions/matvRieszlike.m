function [ nLogL, logLcontr, varargout ] = ...
    matvRieszlike( Sigma_, n, X, varargin )
%MATVRIESZLIKE
%
% USAGE:
%   [ nLogL, score, param ] = matvFlike( Sigma_, n1, n2, X )
%   [ nLogL, score, param ] = matvFlike( [], [], [], X, all_param )
%
% INPUTS:
%   X       - Array (p by p). Symmetric p.d. data matrix.
%   SIGMA_  - Array (p by p). Symmetric p.d. parameter matrix, 
%             regulates covariance. 
%   n1    - Double. First degrees of freedom parameter. 
%   n2    - Double. Second degrees of freedom parameter. 
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
%   [1] Gupta and Nagar (2001) - Matrix Variate Distributions, p.156.
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
% 11.02.2021

[p,~,N] = size(X);
p_ = p*(p+1)/2;
narginchk(3,4);
nargoutchk(0,5);
%% Param
if nargin == 4
    if ~(isempty(Sigma_) && isempty(n))
        error('Cannot input all_param and any of the parameters individually!')
    end
    all_param = varargin{1};
    Sigma_ = ivechchol(all_param(1 : p_));
    n = all_param(p_ + 1 : p_ + p);
end

% Optional Parameter Output
if nargout == 5
    param.Sigma_ = Sigma_;
    param.chol_Sigma_ = vechchol(Sigma_);
    param.n = n;
    param.all = [param.chol_Sigma_; n];
    
    varargout{3} = param;
end
%% Input Checking
chol(Sigma_,'lower'); % Checking if Sigma_ is symmetric p.d.
%% Log-likelihood computation
logLcontr = NaN(N,1);

term1 = -sum(n)/2*log(2);
S = -lgmvgammaln(n./2);
term3 = -loglpwdet(Sigma_,n./2);

log_normalizing_constant = term1 + S + term3;

for ii = 1:N
    term4 = loglpwdet(X(:,:,ii),(n-p-1)./2);
    term5 = -trace(Sigma_\X(:,:,ii))./2;
    
    log_kernel = term4 + term5;

    logLcontr(ii) = log_normalizing_constant + log_kernel;
end
nLogL = -sum(logLcontr);
%% Score computation
if nargout >= 3
    
    score.Sigma_ = NaN(N,p_);
    score.n = NaN(N,p);
    
    for ii = 1:N
        
        % General matrix derivative (ignoring symmetry of Sigma_):
        S = dloglpwdet_dX(Sigma_,-n/2) ...
            + .5*(Sigma_\X(:,:,ii)/Sigma_);
        
        % Accounting for symmetry of Sigma_:
        S = 2*S - diag(diag(S));
        
        % Numerical inaccuracies make S not exactly symmetric
        S = (S+S')./2;
                
        score.Sigma_(ii,:) = vech(S);

        % Score nu
        score.n(ii,:) = NaN(p,1);

    end
    
    varargout{1} = score;

end
%% Fisher Info
if nargout >= 4
    
    fisherinfo.Sigma_ = NaN;
    fisherinfo.n = NaN;
    
    varargout{2} = fisherinfo;
end

end
