function [ nLogL, logLcontr, param, varargout ] = matvFGammamixlike( Sigma_, df_1, df_2, lambda_, data_mat, varargin )
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
%
% Michael Stollenwerk
% michael.stollenwerk@live.com
% 12.07.2020
%
% DEPENDENCIES:
%
%%
[p,~,N] = size(data_mat);
p_ = p*(p+1)/2;
%% Parameters
if nargin == 6
    all_param = varargin{1};
    Sigma_ = ivechchol(all_param(1:p_));
    df_1 = all_param(p_ + 1);
    df_2 = all_param(p_ + 2);
    lambda_ = all_param(p_ + 3);
end

param.Sigma_ = Sigma_;
param.df_1 = df_1;
param.df_2 = df_2;
param.lambda_ = lambda_;
param.all = [vech(chol(Sigma_,'lower'))', df_1, df_2, lambda_];
%% Log-Likelihood 
% logLcontr = NaN(N,1);
% 
% % normalizing constant
% log_norm_constant = ...
%     lambda_ * log(lambda_) ...
%     - gammaln(lambda_) ...
%     - mvbetaln( df_1/2, (df_2 + p - 1)/2, p) ...
%     + (df_2 + p - 1)/2*log(det(Sigma_));
%    
% % kernel
% for ii = 1:N
%     
%     A = data_mat(:,:,ii);
%     
%     log_term1 = (df_1 - p - 1)/2*log(det(A));
%     
%     fun = @(w) ...
%         w.^(lambda_ - 1 - p.^2 - p.*(df_1 - p  - 1)/2) ...
%         .* exp(-lambda_.*w) ...
%         .* arrayfun(@(x) det( Sigma_ + A./x ), w).^(-(df_1 + df_2 + p - 1)./2);
%     term2 = integral( fun, 0, inf );
%     
%     logLcontr(ii) = log_norm_constant + log_term1 + log(term2);
%     
% end

Sigma_rescaled = Sigma_*df_1/(df_2 - 2);
df_2_rescaled = df_2 + p - 1;

% normalizing constant
log_norm_constant = ...
    lambda_ * log(lambda_) ...
    - gammaln( lambda_ ) ...
    - mvbetaln(df_1/2, df_2_rescaled/2, p) ...
    + p * df_1 / 2 * log( df_1 / (df_2_rescaled - p - 1) ) ...
    + (df_1 + 2*df_2_rescaled) / 2 * log(det(Sigma_rescaled));
    
% kernel
logLcontr = NaN(N,1);
% gamma_rnd = gamrnd( lambda_, 1/lambda_, 1e4, 1 );
for ii = 1:N
    
    A = data_mat(:,:,ii);
    
    log_term1 = -(df_2_rescaled + p + 1) / 2 * log(det(A));
    
    Sig_iA_Sig = Sigma_rescaled/A*Sigma_rescaled;
    
% This is using the Expected value formula. Comment in gamma_rnd above and 
% comment out the first two terms of log_norm_constant above!
%     X = gamma_rnd.^(-p*df_1/2) ...
%         .* arrayfun(...
%             @(x) det( Sig_iA_Sig + df_1 / (df_2_star - p - 1) ./ x .* Sigma_star ), ...
%             gamma_rnd ...
%         ).^( -(df_1 + df_2_star) / 2 );
%     
%     term2 = mean(X);

%   Here I tried to make the integral more stable.
%     fun = @(w) ...
%         exp( ...
%             ( lambda_ - 1 - p * df_1 / 2 ).*log(w) ...
%             - lambda_ .* w  ...
%             -(df_1 + df_2_star) / 2 ...
%             .* arrayfun(...
%                 @(x) log(det( Sig_iA_Sig + df_1 / (df_2_star - p - 1) ./ x .* Sigma_star )), ...
%                 w ...
%             ) ...
%         );
%     term2 = integral( fun, 0, inf );
    
%   Here I tried to put all variables under the integral (and make it more stable).
%     fun = @(w) ...
%         exp(log_norm_constant) * exp(log_term1) * ...
%         w.^( lambda_ - 1 - p * df_1 / 2 ) ...
%         .* exp( -lambda_ .* w ) ...
%         .* arrayfun(...
%             @(x) det( Sig_iA_Sig + df_1 / (df_2_star - p - 1) ./ x .* Sigma_star ), ...
%             w ...
%         ).^( -(df_1 + df_2_star) / 2 );
%     logLcontr(ii) = log(integral( fun, 0, inf ));
    
%     term2 = integral( fun, 0, inf );    
    
%   Here I programmed the function fun at the end of the code.
%   term2 = integral( @(x) fun(x,Sig_iA_Sig), 0, 1e3/lambda_, 'RelTol', 0, 'AbsTol', 1e-12 );
    
    fun = @(w) ...
        w.^( lambda_ - 1 - p * df_1 / 2 ) ...
        .* exp( -lambda_ .* w ) ...
        .* arrayfun(...
            @(x) det( Sig_iA_Sig + df_1 / (df_2_rescaled - p - 1) ./ x .* Sigma_rescaled ), ...
            w ...
        ).^( -(df_1 + df_2_rescaled) / 2 );
	term2 = integral( fun, 0, inf, 'RelTol', 0, 'AbsTol', 1e-12 );
    
    logLcontr(ii) = log_norm_constant + log_term1 + log(term2);
        
end
nLogL = -sum(logLcontr);

%% Score
if nargout == 4
    score.Sigma_ = NaN(p);
    score.df_1 = NaN;
    score.df_2 = NaN;
    score.lambda = NaN;
    varargout = score;
end

% %% Nested fun
% function y = fun(w,Sig_iA_Sig) 
%     term1_ = w.^( lambda_ - 1 - p * df_1 / 2 );
%     term2_ = exp( -lambda_ .* w );
%     term3_ = NaN(1,length(w));
%     for jj = 1:length(w)
%         term3_(jj) = ...
%            det( Sig_iA_Sig + df_1 / (df_2_star - p - 1) ./ w(jj) .* Sigma_star );
%     end
%     y = term1_.*term2_.*(term3_.^( -(df_1 + df_2_star) / 2 ));
% end

end
