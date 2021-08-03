clear
clc
%%
T = 4e3;

Omega_ = repmat([.35,.5;.5,.75],1,1,T);
Sigma_ = repmat(eye(2,2),1,1,T);
df_F_1 = ones(T,1)*10;
df_F_2 = ones(T,1)*15;

param = [df_F_1(1); df_F_2(1); vech(chol(Sigma_(:,:,1),'lower')); vech(chol(Omega_(:,:,1),'lower'))]';

options = optimoptions('fminunc','Display','final','MaxFunEval',1e5);
eparam = NaN(20,8);
parfor ii = 1:20
    data_mat = mncFrnd( df_F_1, df_F_2, Sigma_, Omega_);

    fun_ = @(x) fun(x,data_mat);

    eparam(ii,:) = fminunc(fun_,[5,10,1,0,1,1,0,1],options);
end

eparam

for ii = 1:20
    eSigma(:,:,ii) = ivech(eparam(ii,3:5),'lower')*ivech(eparam(ii,3:5),'lower')';
    eOmega(:,:,ii) = ivech(eparam(ii,6:8),'lower')*ivech(eparam(ii,6:8),'lower')';
    eEV(:,:,ii) = eSigma(:,:,ii) + eOmega(:,:,ii)/eparam(ii,1);
end
mean(eSigma,3)
mean(eOmega,3)
mean(eEV,3)
%%
function y = fun(x,data_mat)
if size(x,1) < size(x,2)
    x = x';
end    
T = size(data_mat,3);
df_F_1 = ones(T,1)*x(1);
df_F_2 = ones(T,1)*x(2);
Sigma_ = ivech(x(3:5),'lower');
Sigma_ = repmat(Sigma_*Sigma_',1,1,T);
Omega_ = ivech(x(6:8),'lower');
Omega_ = repmat(Omega_*Omega_',1,1,T);
y = mncFlike(data_mat, df_F_1, df_F_2, Sigma_, Omega_, 100);
end
