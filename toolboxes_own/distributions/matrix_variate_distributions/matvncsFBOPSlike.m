function [ nLogL, logLcontr, varargout ] = ...
    matvncsFBOPSlike( Sigma_, m, df_1, df_2, X, varargin )
%MATVNCSFBOPSLIKE Negative log-likelihood and score for the matrix-variate 
% non-central spiked F distr.
%
% USAGE:
%  [ nLogL, score, param ] = matvncFBOPSlike( Sigma_, Omega_, df_1, df_2, mhg_precision, X)
%  [ nLogL, score, param ] = matvncFBOPSlike( [], [], [], [], mhg_precision, X, all_param)
%
% INPUTS:
%   SIGMA_        - Array (p by p). Symmetric p.d. parameter matrix.
%   OMEGA_        - Array (p by p). Symmetric p.d. non-centrality matrix.
%   DF_1          - Double. First degrees of freedom parameter. 
%   DF_2          - Double. Second degrees of freedom parameter. 
%   MHG_PRECISION - Number of Jack Functions used to approximate the matrix
%                   valued hypergeometric function. This input is inversly
%                   related to computing time.
%   X             - Array (p by p). Symmetric p.d. data matrix.
%   VARARGIN      - Array (p*(p+1) + 2 by 1). All parameters in one array.
%                   [vech(chol(SIGMA_,'lower'));vech(chol(OMEGA_,'lower'));DF_1;DF_2]
%
% OUTPUTS:
%   NLOGL        - Negative log-likelihood value
%   SCORE        - Struct. Fields as parameter names. Contains derivatives
%                  of log-likelihood w.r.t. parameters.
%   VARARGOUT    - Struct. Optional. Parameters.
%
% COMMENTS:
%
% REFERENCES:
%      [1] Bodnar, Okhrin, Parolya and Stollenwerk (2020)
%
% DEPENDENCIES:
%   MVBETALN IVECHCHOL
%
% See also MVBETALN IVECHCHOL
%
% Michael Stollenwerk
% michael.stollenwerk@live.com
% 16.06.2020

[p,~,N] = size(X);
p_ = p*(p+1)/2;
narginchk(5,6);
nargoutchk(0,5);
%% Param
if nargin == 6
    if ~(isempty(Sigma_) && isempty(m) && isempty(df_1) && isempty(df_2))
        error('Cannot input all_param and any of the parameters individually!')
    end
    all_param = varargin{1};
    if size(all_param,1) < size(all_param,2)
        all_param = all_param';
    end
    Sigma_ = ivechchol(all_param(1 : p_));
    m = all_param(p_ + 1: p_ + p);
    df_1 = all_param(p_ + p + 1);
    df_2 = all_param(p_ + p + 2);
end

% Optional Parameter Output
if nargout == 5
    
    param.Sigma_ = Sigma_;
    param.chol_Sigma_ = vechchol(Sigma_);
    param.m = m;
    param.df_1 = df_1;
    param.df_2 = df_2;
    param.all = [param.chol_Sigma_; param.m; df_1; df_2];

    varargout{3} = param;
    
end
%% Log-likelihood computation
logLcontr = NaN(N,1);

y_hg = NaN(N,1);

invSigma = inv(Sigma_);
invY = NaN(p,p,N);
for ii = 1:N
    
    c1 = df_1 / (df_2 - p - 1);
    
    Y = Sigma_/X(:,:,ii)*Sigma_ + c1*Sigma_;
    invY(:,:,ii) = inv(Y);
    if any(isnan(invY(:))) || any(isinf(invY(:)))

        nLogL = NaN;
        score.Sigma_ = NaN(N,p_);
        score.m = NaN(N,p);
        score.df_1 = NaN(N,1);
        score.df_2 = NaN(N,1);
        varargout{1} = score;

        fisherinfo.Sigma_ = NaN;
        fisherinfo.m = NaN;
        fisherinfo.df_1 = NaN;
        fisherinfo.df_2 = NaN;

        varargout{2} = fisherinfo;
        return
    end     
    
    % matrix-F part
    term1 = -mvbetaln(df_1/2, df_2/2, p);
    term2 = p*df_1/2 * log(c1);
    term3 = -(df_2 + p + 1)/2*log(det(X(:,:,ii)));    
    term4 = (df_1 + 2*df_2)/2*log(det(Sigma_));
    term5 = -(df_1 + df_2)/2*log(det(Y));

    mFpart = term1 + term2 + term3 + term4 + term5;
    logLcontr(ii) = mFpart;

    if any(m(:) ~= 0)
        
        % noncentral part
        term6 = -.5.*(m'*invSigma*m);
        
        y_hg(ii) = ...
            hypergeom11_my( ...
                (df_1 + df_2)/2, df_1/2, .5.*c1.*(m'*invY(:,:,ii)*m),100 ...
            );
        term7 = log(y_hg(ii));

        ncpart = term6 + term7;
        logLcontr(ii) = mFpart + ncpart;
    end    
end
nLogL = -sum(logLcontr);

%% Score computation
if nargout >= 3
    
    D = Dmatrix(p);
    
    scoreSig = NaN(N,p_);
    score_m = NaN(N,p);
    score_df_1 = NaN(N,1);
    score_df_2 = NaN(N,1);    
    for ii = 1:N
        % Score Sigma_   
        S_Sig = ...
            (df_1/2 + df_2)*invSigma ...
            - (df_1 + df_2)/2*(invY(:,:,ii)*Sigma_/X(:,:,ii) + ...
                               X(:,:,ii)\Sigma_*invY(:,:,ii) + ...
                               c1*invY);   
        
        if any(m(:) ~= 0)
            
            const2 = df_2/2*hypergeom11_my( ...
                             (df_1 + df_2)/2+1, df_1/2+1, .5.*c1.*(m'*invY(:,:,ii)*m), 100 ...
                         )/y_hg(ii);
                     
            S_Sig = S_Sig ...
                + .5.*invSigma*(m*m')*invSigma ...
                - const2*c1/2*(invY(:,:,ii)*(m*m')*invY(:,:,ii)*Sigma_/X(:,:,ii) + ...
                                          X(:,:,ii)\Sigma_*invY(:,:,ii)*(m*m')*invY(:,:,ii) + ...
                                          c1*invY(:,:,ii)*(m*m')*invY(:,:,ii));
                                      
            score_m(ii,:) = -invSigma*m + const2*invY*m;
                                      
        end
        
        S_Sig = S_Sig(:); % Apply vec operator.
        scoreSig(ii,:) = D'*S_Sig; % Accounting for symmetry of Sigma_.

    end
    
    score.Sigma_ = scoreSig;
    score.m = score_m;
    score.df_1 = score_df_1;
    score.df_2 = score_df_2;
    
    varargout{1} = score;
end
%% Fisher Info
if nargout >= 4
    
    fisherinfo.Sigma_ = NaN;
    fisherinfo.m = NaN;
    fisherinfo.df_1 = NaN;
    fisherinfo.df_2 = NaN;
    
    varargout{2} = fisherinfo;
end
end
