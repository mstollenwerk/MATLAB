function [ nLogL, logLcontr, Sigma_, Omega_ ] = cancF_likeRec(intrcpt, arch_param, garch_param, war_matrix, df_F_1, df_F_2, data_RC)
%% Input Checking
% Will be added later
[k,~,T] = size(data_RC);
%% Data Storage
Sigma_ = NaN(k,k,T+1);
Omega_ = NaN(k,k,T+1);
 %% Recursion
for tt=1:T+1
    if (tt-1) <= 0
        Sigma_(:,:,tt) = intrcpt + arch_param*mean(data_RC,3) + garch_param*mean(data_RC,3);
        Omega_(:,:,tt) = war_matrix*mean(data_RC,3)*war_matrix';
    else
        Sigma_(:,:,tt) = intrcpt + arch_param*data_RC(:,:,tt-1) + garch_param*Sigma_(:,:,tt-1);
        Omega_(:,:,tt) = war_matrix*data_RC(:,:,tt-1)*war_matrix';
    end
end
%% Likelihood
[ nLogL, logLcontr ] = mncFlike(...
    data_RC, ones(T,1)*df_F_1, ones(T,1)*df_F_2,...
    Sigma_(:,:,1:end-1), Omega_(:,:,1:end-1), 10);
