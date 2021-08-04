function [ nLogL, logLcontr, varargout ] = ...
    matvgigWishlike( ...
        Sigma_, n, lambda_, xi_, psi_, X, varargin ...
    )
%
%
% USAGE:    [ nLogL, logLcontr, varargout ] = ...
%               matvqghypberboliclike( ...
%                   Sigma_, df_n, lambda_, xi_, psi_, data_mat, varargin ...
%               )
%   
%
% INPUTS:
%   
%
% OUTPUTS:
%   
%  See also 
%
% COMMENTS:
%   
% REFERENCES:
%      [1]                    
%
% Michael Stollenwerk
% michael.stollenwerk@live.com
% 25.09.2020
%
% DEPENDENCIES:
%

[p,~,N] = size(X);
p_ = p*(p+1)/2;
narginchk(6,7);
nargoutchk(0,5);
%% Parameters
if nargin == 7
    all_param = varargin{1};
    Sigma_ = ivechchol(all_param(1:p_));
    n = all_param(p_ + 1);
    lambda_ = all_param( p_ + 2);
    xi_ = all_param( p_ + 3);
    psi_ = all_param( p_ + 4);
end

if nargout == 5
    param.Sigma_ = Sigma_;
    param.df_n = n;
    param.lambda_ = lambda_;
    param.xi_ = xi_;
    param.psi_ = psi_;
    param.all = [vech(chol(Sigma_,'lower'))', n, lambda_, xi_, psi_]';
    
    varargout{3} = param;
end
%% Input Checking
cholSig = chol(Sigma_,'lower'); % Checking if Sigma_ is symmetric p.d.
%% Log-Likelihood
logLcontr = NaN(N,1);

log_norm_const = ...
    (p*n)/2*log(pi) ...
    - mvgammaln( n/2, p ) ...
    + (p*n)/2 * log(psi_/2/pi) ...
    - lambda_/2 * log(xi_*psi_) ...
    - n/2*logdet(Sigma_);

term4 = - log(besselk( sym(lambda_), sym(sqrt(xi_*psi_)) ));
term4 = double(term4);

log_norm_const = log_norm_const + term4;

for ii = 1:N
    
    A = X(:,:,ii);
    q = xi_ + trace(Sigma_\A);
    
    term1 = (n - p - 1)/2*log(det(A));
    term2 = log(besselk( sym(lambda_ - (p*n)/2), sym(sqrt(psi_*q)) ));
    term2 = double(term2);
    term3 = - ((p*n)/2 - lambda_) * log(sqrt(psi_*q));
    log_kernel = term1 + term2 + term3;
    
    logLcontr(ii) = log_norm_const + log_kernel;

end
nLogL = -sum(logLcontr);
%% Scores
if nargout >= 3
    invSig = inv(Sigma_);

    score.Sigma_ = NaN(N,p_);
    score.n = NaN(N,1);
    score.lambda = NaN(N,1);
    score.xi = NaN(N,1);
    score.psi = NaN(N,1);
    
    for ii = 1:N
        qq = invSig*X(:,:,ii)*invSig;
        S = ...
            n*invSig ...
            + (   besselk( lambda_ - (p*n)/2 - 1, sqrt(psi_*q) ) ...
                + besselk( lambda_ - (p*n)/2 + 1, sqrt(psi_*q) )    ) ...
              /exp(term2)/sqrt(psi_*q)*psi_*qq ...
            + (p*n/2 - lambda_)/q*qq;
        S = .5.*S;

        % Accounting for symmetry of Sigma_:
        S = 2*S - diag(diag(S));
        
        % Numerical inaccuracies make S not exactly symmetric
        S = (S+S')./2;
                
        score.Sigma_(ii,:) = NaN;%vech(S);

        % Score dfs
        score.n(ii) = NaN;
        score.lambda(ii) = NaN;
        score.xi(ii) = NaN;
        score.psi(ii) = NaN;

    end
    
    varargout{1} = score;

end
%% Fisher Info
if nargout == 4
    
    varargout{2} = NaN;
    
end
end
