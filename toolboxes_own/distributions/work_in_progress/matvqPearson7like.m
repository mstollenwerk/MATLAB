function [ nLogL, logLcontr, param, varargout ] = matvqPearson7like( Sigma_, n, q, r, X, varargin )
%BESSELWISHLIKE
%
% USAGE:
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
%      [1] Based on Bodnar, Gupta and Varga p.54
%
% Michael Stollenwerk
% michael.stollenwerk@live.com
% 12.07.2020
%
% DEPENDENCIES:
%
%%
[p,~,N] = size(X);
p_ = p*(p+1)/2;
narginchk(5,6); %%%%%%%
nargoutchk(1,5);
%% Parameters
if nargin == 6
    if ~(isempty(Sigma_) && isempty(n) && isempty(q) && isempty(r))
        error('Cannot input all_param and any of the parameters individually!')
    end    
    all_param = varargin{1};
    Sigma_ = ivechchol(all_param(1:p_));
    n = all_param(p_ + 1);
    q = all_param(p_ + 2);
    r = all_param(p_ + 3);
end

% Optional Paramter Output
if nargout == 5
    param.Sigma_ = Sigma_;
    param.n = n;
    param.q = q;
    param.r = r;
    param.all = [vech(chol(Sigma_,'lower')); n; q; r];
    
    varargout{3} = param;    
end
%% Log-Likelihood 
% normalizing constant
log_norm_constant = ...
      gammaln( q ) ...
    - gammaln( q - p*n/2 ) ...
    - p*n/2*log( r ) ...
    - n/2*log(det( Sigma_ )) ...
    - mvgammaln( n/2, p );

% kernel
log_kernel = NaN(N,1);
for ii = 1:N
    
    A = X(:,:,ii);
      
    log_kernel(ii) = (n - p - 1)/2*log(det( A )) ...
        - q*log( 1 + trace(Sigma_\A)/r );
    
end
logLcontr = log_norm_constant + log_kernel;
nLogL = -sum(logLcontr);
%% Score computation
if nargout >= 3
    
    scoreSig = NaN(N,p_);
    score_n = NaN(N,1);
    score_q = NaN(N,1);
    score_r = NaN(N,1);
    
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
    
    score.Sigma_ = scoreSig;
    score.n = score_n;
    score.q = score_q;
    score.r = score_r;
    
    varargout{1} = score;

end
%% Fisher Info
if nargout >= 4
    
    fisherinfo.Sigma_ = NaN(p_,p_,N);
    fisherinfo.n = NaN(N,1);
    fisherinfo.q = NaN(N,1);
    fisherinfo.r = NaN(N,1);
    
    varargout{2} = fisherinfo;
end

end
