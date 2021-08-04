function [ nLogL, logLcontr, varargout ] = matvqtGammaprodlike( Sigma_, df_n, df_t, lambda_, data_mat, varargin )
%
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
%      [1]                    
%
% Michael Stollenwerk
% michael.stollenwerk@live.com
% 14.07.2020
%
% DEPENDENCIES:
%

narginchk(5,6);
nargoutchk(1,5);

[p,~,N] = size(data_mat);
p_ = p*(p+1)/2;
%% Parameters
if nargin == 6
    all_param = varargin{1};
    Sigma_ = ivechchol(all_param(1:p_));
    df_n = all_param(p_ + 1);
    df_t = all_param(p_ + 2);
    lambda_ = all_param(p_ + 3);
end

if nargout == 5
    param.Sigma_ = Sigma_;
    param.df_n = df_n;
    param.df_t = df_t;
    param.lambda_ = lambda_;
    param.all = [vech(chol(Sigma_,'lower'))', df_n, df_t, lambda_];
    
    varargout{3} = param;
end
%% Log-Likelihood
logLcontr = NaN(N,1);

log_norm_const = ...
    gammaln( ( p*df_n + df_t )/2 ) ...
    - mvgammaln(df_n/2,p) ...
    - betaln( lambda_, df_t/2 ) ...
    + p*df_n/2*log( lambda_ / df_t ) ...
    - df_n/2*log(det(Sigma_));
    
for ii = 1:N
    
    A = data_mat(:,:,ii);

%     try
%     log_kern = ...
%         (df_n - p - 1) / 2 * log(det(A)) ...
%         + log( ...
%             kummerU_taylorb( ...
%                 ( p*df_n + df_t )/2, ...
%                 p*df_n/2 - lambda_ + 1, ...
%                 lambda_ / df_t * trace(Sigma_\A), ...
%                 eps ...
%             ) ...            
%         );
%     catch
%     log_kern = ...
%         (df_n - p - 1) / 2 * log(det(A)) ...
%         + log( ...
%             kummerU( ...
%                 ( p*df_n + df_t )/2, ...
%                 p*df_n/2 - lambda_ + 1, ...
%                 lambda_ / df_t * trace(Sigma_\A) ...
%             ) ...            
%         );
%     end
    
    % This is my c implemented kummerU function based on the arb.c library.
    log_kernel = ...
        (df_n - p - 1) / 2 * log(det(A)) ...
        + log( ...
            kummerU_my( ...
                ( p*df_n + df_t )/2, ...
                p*df_n/2 - lambda_ + 1, ...
                lambda_ / df_t * trace(Sigma_\A), ...
                1000 ...
            ) ...            
        );

    logLcontr(ii) = log_norm_const + log_kernel;
    
%     if isinf(logLcontr(ii))
%         nLogL = logLcontr(ii);
%         return
%     end
%     
%     if isnan(logLcontr(ii))
%         nLogL = logLcontr(ii);
%         return
%     end
    
end
nLogL = -sum(logLcontr);
%% Scores
if nargout >= 3
    invSig = inv(Sigma_);
    
    scoreSig = NaN(N,p_);
    score_df_n = NaN(N,1);
    score_df_t = NaN(N,1);
    for ii = 1:N
        
        A = data_mat(:,:,ii);
       
        Nenner = kummerU_my( ...
                ( p*df_n + df_t )/2, ...
                p*df_n/2 - lambda_ + 1, ...
                lambda_ / df_t * trace(Sigma_\A), ...
                1000 ...
            );
        Zaehler = kummerU_my( ...
                ( p*df_n + df_t )/2 + 1, ...
                p*df_n/2 - lambda_ + 1 + 1, ...
                lambda_ / df_t * trace(Sigma_\A), ...
                1000 ...
            );
        S = - df_n/2*invSig ...
            + (df_t + p*df_n)/ 2 * Zaehler/Nenner * (Sigma_\A/Sigma_);
        
        % Apply vec operator.
        S = S(:);
        
        % Accounting for symmetry of Sigma_:
        D = Dmatrix(p);
        scoreSig(ii,:) = D'*S;

        % The score below are also easy to get quering wolframalpha.com 
        % with eg "d/da (log(Gamma(1/2 (a+ p n))))".
        % I am just too lazy to write them down right now.
        score_df_n(ii) = NaN;
        score_df_t(ii) = NaN;
    
    end
    
    score.Sigma_ = scoreSig;
    score.df_n = score_df_n;
    score.df_t = score_df_t;
    
    varargout{1} = score;
    
end
%% Fisher Info
if nargout == 4
        
    % THIS IS WRONG. IT IS THE SCALING OF THE WISHART.
    fisherInfo.Sigma_ = df_n/2*D'*kron(invSig,invSig)*D;
    
    varargout{2} = fisherInfo;
end
end
