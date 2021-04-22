function [ data_RC, Sigma_, Omega_ ] = cancF_rnd( param, N )
error('Not done yet. This is dummy code.')
df_F_1 = 10;
df_F_2 = 15;
intrcpt = param.intrcpt;
arch_param = param.arch_param;
garch_param = param.garch_param;
war_matrix = param.war_matrix;
%% Input Checking
% Furter input checks will be added later
[K,~,~] = size(intrcpt);
%% Initialization
burn_in = 1e3;
N__ = N + burn_in;

Sigma_ = NaN(K,K,N__);
Sigma_(:,:,1) = intrcpt;

Omega_ = NaN(K,K,N__);
Omega_(:,:,1) = war_matrix*intrcpt;

data_RC = NaN(K,K,N__);
data_RC(:,:,1) = mncsFrnd(df_F_1, df_F_2, Sigma_(:,:,1), Omega_(:,:,1));
%% Recursion
for tt=2:N__
    Sigma_(:,:,tt) = intrcpt + arch_param*data_RC(:,:,tt-1) + ...
        garch_param*Sigma_(:,:,tt-1);
    Omega_(:,:,tt) = war_matrix*data_RC(:,:,tt-1)*war_matrix';
    
    data_RC(:,:,tt) = mncsFrnd(df_F_1, df_F_2, Sigma_(:,:,tt), Omega_(:,:,tt));
end
%% Deduct Burn In
Sigma_ = Sigma_(:,:,burn_in+1:end);
Omega_ = Omega_(:,:,burn_in+1:end);
data_RC = data_RC(:,:,burn_in+1:end);
