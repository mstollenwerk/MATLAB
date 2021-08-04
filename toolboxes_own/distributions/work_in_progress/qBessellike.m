function [ nLogL, logLcontr ] = qBessellike( data_mat, df_1, df_2, df_3, Sigma_ )
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
%      [1] Based on Bodnar, Gupta and Varga but apparent mistake
%      (tr(Q)^{q/2}) corrected.
%
% Michael Stollenwerk
% michael.stollenwerk@live.com
% 18.02.2020
%
% DEPENDENCIES:
%
[p,~,N] = size(data_mat);
%% Log-likelihood computation
logLcontr=NaN(N,1);

for ii = 1:N
    
    trace_invSig_dta = trace(Sigma_(:,:,ii)\data_mat(:,:,ii));
    
    term1 = -(df_2(ii)+p*df_1(ii)-1)*log(2);
    term2 = -(p*df_1(ii)+df_2(ii))*log(df_3(ii));
    term3 = -mvgammaln(df_1(ii)/2,p) - gammaln(df_2(ii) + p*df_1(ii)/2);
    term4 = -df_1(ii)/2*log(det(Sigma_(:,:,ii)));
    log_kernel = df_2(ii)/2*log(trace_invSig_dta) + ...
                 (df_1(ii)-p-1)/2*log(det(data_mat(:,:,ii))) + ...
                 log(besselk(df_2(ii),sqrt(trace_invSig_dta)/df_3(ii)));
             
    logLcontr(ii) = term1 + term2 + term3 + term4 + log_kernel;    
end

nLogL = -sum(logLcontr);
end
