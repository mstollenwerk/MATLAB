function [ nLogL, logLcontr, varargout ] = ...
    matvncFBOPSlike( Sigma_, Omega_, df_1, df_2, mhg_precision, X, varargin )
%MATVNCFBOPSLIKE Negative log-likelihood and score for the matrix-variate 
% non-central F distr.
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
%   Description of MHG:
%       Code by Koev and Edelman (2006)
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
%   MVBETALN MHG.C MHG IVECHCHOL
%
% See also MVBETALN MHG.C MHG IVECHCHOL
%
% Michael Stollenwerk
% michael.stollenwerk@live.com
% 16.06.2020

[p,~,N] = size(X);
p_ = p*(p+1)/2;
narginchk(6,7);
nargoutchk(1,5);
%% Param
if nargin == 7
    if ~(isempty(Sigma_) && isempty(Omega_) && isempty(df_1) && isempty(df_2))
        error('Cannot input all_param and any of the parameters individually!')
    end
    all_param = varargin{1};
    if size(all_param,1) < size(all_param,2)
        all_param = all_param';
    end
    Sigma_ = ivechchol(all_param(1 : p_));
    Omega_ = ivechchol(all_param(p_ + 1: p_ + p_));
    df_1 = all_param(p_ + p_ + 1);
    df_2 = all_param(p_ + p_ + 2);
end

% Optional Parameter Output
if nargout == 5
    
    param.Sigma_ = Sigma_;
    param.chol_Sigma_ = vechchol(Sigma_);
    param.Omega_ = Omega_;
    param.chol_Omega_ = vechchol(Omega_);
    param.df_1 = df_1;
    param.df_2 = df_2;
    param.all = [param.chol_Sigma_; param.chol_Omega_; df_1; df_2];

    varargout{3} = param;
    
end
%% Boundary Cases
% Omega_ == 0 is built into Log-likelihood computation below.
%% Input Checking
% if X ~= X'
%     error('data_mat is not quadratic symmetric.')
% end
% 
% if Sigma_ ~= Sigma_'
%     error('Sigma_ is not quadratic symmetric.')
% end
% 
% if Omega_ ~= Omega_'
%     error('Omega_ is not quadratic symmetric.')
% end
% 
% if any(size(Sigma_) ~= size(X)) || any(size(Omega_) ~= size(X))
%     error('Size of Sigma, Omega_ and data_mat does not match')
% end

if isempty(mhg_precision)
    mhg_precision = 50;
end

% disp('MHG precision parameter (number of jack functions to be)')
% disp(['calculated is set to ', num2str(mhg_precision), '. This is strongly related to'])
% disp('computing time.')

%% Log-likelihood computation
logLcontr = NaN(N,1);

y_mhg = NaN(N,1);
for ii = 1:N
    
    c1 = df_1 / (df_2 - p - 1);
    
    Y = Sigma_/X(:,:,ii)*Sigma_ + c1*Sigma_;
    
    % matrix-F part
    term1 = -mvbetaln(df_1/2, df_2/2, p);
    term2 = p*df_1/2 * log(c1);
    term3 = -(df_2 + p + 1)/2*log(det(X(:,:,ii)));    
    term4 = (df_1 + 2*df_2)/2*log(det(Sigma_));
    term5 = -(df_1 + df_2)/2*log(det(Y));

    mFpart = term1 + term2 + term3 + term4 + term5;
    logLcontr(ii) = mFpart;

    if any(Omega_(:) ~= 0)
        
        % noncentral part
        term6 = trace(-.5.*(Sigma_\Omega_));
        % Only the nonzero eigenvalues are taken (if not all -> spiked ncF).
        x_mhg = eig(.5.*c1.*(Omega_/Y));
        if all(x_mhg == 0)
            x_mhg = 0;
        else
            x_mhg = x_mhg(x_mhg ~= 0);
        end
        try
            y_mhg(ii) = mhg(...
                mhg_precision, 2, ...
                (df_1 + df_2)/2, df_1/2, x_mhg ...
            );
            term7 = log(y_mhg(ii));
        catch
            
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

        ncpart = term6 + term7;
        logLcontr(ii) = mFpart + ncpart;
    end    
end
nLogL = -sum(logLcontr);

%% Score computation
if nargout >= 3
    
    invSigma = inv(Sigma_);

    scoreSig = NaN(N,p_);
    scoreOm = NaN(N,p_);
    score_df_1 = NaN(N,1);
    score_df_2 = NaN(N,1);    
    for ii = 1:N
    % Score Sigma_

    % Nestors letzte falsche Formel:
    %     c2 = (df_1 + df_2) / (df_2 - p - 1);
    %     X_SigOm = invX*Sigma_*Omega_;
    %     
    %     Z = - df_1*invSigma ...
    %          + df_1*c2*invY ...
    %          + invSigma*Omega_*invSigma ...
    %          - c2*(X_SigOm + X_SigOm' + c1*Omega_);  
    %     Z = .5.*(Z+Z'); % Symmetric Sigma and Omega in -> Symmetric Score out.
    %     Z = .5*Z;

    % Oben Wishart approx unten exakt oder, wenn Q = mean_ncW, dann alte
    % falsche Formel, wenn Q = (c1 + c2)*Sigma_ncW + Omega_ncW, dann andere
    % alte falsche Formel:
    
        Y = Sigma_/X(:,:,ii)*Sigma_ + c1*Sigma_;
        invY = inv(Y);
        if any(isnan(Y(:))) || any(isinf(Y(:))) || any(isnan(invY(:))) || any(isinf(invY(:)))

            nLogL = NaN;
            score.Sigma_ = NaN(N,p_);
            score.Omega_ = NaN(N,p_);
            score.nu1 = NaN(N,1);
            score.nu2 = NaN(N,1);
            varargout{1} = score;
            
            fisherinfo.Sigma_ = NaN;
            fisherinfo.Omega_ = NaN;
            fisherinfo.df_1 = NaN;
            fisherinfo.df_2 = NaN;

            varargout{2} = fisherinfo;
            return
        end    
        if any(Omega_(:) ~= 0)

            Sigma_ncW = inv(Y);
            Omega_ncW = c1*(Y\Omega_/Y);
            Theta_ncW = Sigma_ncW\Omega_ncW;
            mean_ncW = df_1.*Sigma_ncW + Omega_ncW;

            % Use central wishart Approximation in Nominator and Formula 3.5.19
            % Gupta Nagar in Denominator.
            Q = mean_ncW.*(1 + df_2/df_1).*det(eye(p) + Theta_ncW./df_1)^(df_2/2) ...
                ./exp(trace(-.5*Theta_ncW))./y_mhg(ii);
        %     Q = mean_ncW;
        %     Q = (df_1 + df_2)*Sigma_ncW + Omega_ncW;
            last_term = -.5*( Q*Sigma_/X(:,:,ii) + X(:,:,ii)\Sigma_*Q + c1*Q );

            S_Sig = (.5* df_1 + df_2)*invSigma ...
                + .5* invSigma*Omega_*invSigma ...
                + last_term;
        else
            Sigma_ncW = inv(Y);
            Q = Sigma_ncW.*(df_1 + df_2);
            last_term = -.5*( Q*Sigma_/X(:,:,ii) + X(:,:,ii)\Sigma_*Q + c1*Q );

            S_Sig = (.5* df_1 + df_2)*invSigma + last_term;
        end

    % % This is my try using the proportionality of E(detVh*V) to E(V).   
    % if any(Omega_(:) ~= 0)
    % 
    %     inv_Sigma_ncW = Sigma_ / X * Sigma_ + c1 * Sigma_;
    %     Sigma_ncW = inv(inv_Sigma_ncW);
    %     Omega_ncW = c1*(inv_Sigma_ncW\Omega_)/inv_Sigma_ncW;
    %     Theta_ncW = Sigma_ncW\Omega_ncW;
    %     mean_ncW = df_1.*Sigma_ncW + Omega_ncW;
    %     
    %     term1 = (.5* df_1 + df_2)*invSigma;
    %     term2 = .5* invSigma*Omega_*invSigma;
    %     
    %     term3_c = -.5*( mean_ncW*Sigma_*invX + invX*Sigma_*mean_ncW + c1*mean_ncW );
    %     
    %     mcrep = 1e4;
    %     A = matvncFBOPSrnd(df_1, df_2, Sigma_, Omega_, mcrep);
    %     c = NaN(mcrep,1);
    %     for ii = 1:mcrep
    %         invA = inv(A(:,:,ii));
    %         
    %         term3_c_rnd = ...
    %             -.5*( mean_ncW*Sigma_*invA + invA*Sigma_*mean_ncW + c1*mean_ncW );
    %         
    %         c(ii) = -mean(mean((term1 + term2)./term3_c_rnd));
    %     end
    % 
    %     Z = term1 + term2 + mean(c).*term3_c;
    %     
    % else
    %     inv_Sigma_ncW = Sigma_ / X * Sigma_ + c1 * Sigma_;
    %     Sigma_ncW = inv(inv_Sigma_ncW);
    %     Q = Sigma_ncW.*(df_1 + df_2);
    %     term3 = -.5*( Q*Sigma_*invX + invX*Sigma_*Q + c1*Q );
    % 
    %     Z = (.5* df_1 + df_2)*invSigma + term3;
    % end

    % Exakte Formel aber Zähler durch noncentral Wishart Monte Carlo
    % approximert:
    %     matvncWish_rnd = matvncWishrnd( ...
    %         df_1, ...
    %         Sigma_ncW, ...
    %         Theta_ncW, ...
    %         1e6 ...
    %     );
    %     mcrep = size(normrnd_,2);
    % 
    %     nom = NaN(p,p,mcrep);
    %     for jj = 1:mcrep
    % 
    %         V = matvncWish_rnd(:,:,jj);
    % 
    %         Z_V = -.5*( V*Sigma_*invRC + invRC*Sigma_*V + c*V );
    %         detV = det(V);
    % 
    %         nom(:,:,jj) = Z_V * detV^(df_2/2);
    % 
    %     end
    % 
    %     avg_Z_VdetVh = mean(nom,3);
    %     Z = (.5* df_1 + df_2)*invSigma ...
    %         + .5* invSigma*Omega_*invSigma ...
    %         + avg_Z_VdetVh./mean_detVh;

        % Apply vec operator.
        S_Sig = S_Sig(:);
        
        % Accounting for symmetry of Sigma_:
        D = Dmatrix(p);
        scoreSig(ii,:) = D'*S_Sig;

    % Score Omega
        if any(Omega_(:) ~= 0)

            Omega_ncW = c1*(Omega_/Y*Omega_);
            Theta_ncW = Omega_\Omega_ncW;
            mean_ncW = df_1.*Omega_ + Omega_ncW;
            invOmega = inv(Omega_);

            % Use central wishart Approximation in Nominator and Formula 3.5.19
            % Gupta Nagar in Denominator.
            x_mhg = eig(.5.*(Omega_\Omega_ncW));
            x_mhg = x_mhg(x_mhg ~= 0);
            last_term = ...
                mean_ncW.*(1 + df_2/df_1).*det(eye(p) + Theta_ncW./df_1)^(df_2/2) ...
                ./exp(trace(-.5*Theta_ncW))...
                ./mhg( mhg_precision, 2, (df_1 + df_2)/2, df_1/2, x_mhg );

            S_Om = - .5* (df_1 + df_2)*invOmega ...
                      - .5* invSigma ...
                      + .5* invOmega*last_term*invOmega;

            % Apply vec operator.
            S_Om = S_Om(:);

            % Accounting for symmetry of Sigma_:
            D = Dmatrix(p);
            scoreOm(ii,:) = D'*S_Om;  
        else
            scoreOm(ii,:) = zeros(p_,1);
        end   

        % Score nu
        score_df_1(ii) = NaN;
        score_df_2(ii) = NaN; 
    end
    
    score.Sigma_ = scoreSig;
    score.Omega_ = scoreOm;
    score.df_1 = score_df_1;
    score.df_2 = score_df_2;
    
    varargout{1} = score;
end
%% Fisher Info
if nargout >= 4
    
    fisherinfo.Sigma_ = NaN;
    fisherinfo.Omega_ = NaN;
    fisherinfo.df_1 = NaN;
    fisherinfo.df_2 = NaN;
    
    varargout{2} = fisherinfo;
end
end
