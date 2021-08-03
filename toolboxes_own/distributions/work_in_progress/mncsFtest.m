clear
clc
%%
T = 2e3;

Omega_ = repmat([.5,.75]'*[.5,.75],1,1,T);
Sigma_ = randn(2);
Sigma_ = repmat(Sigma_*Sigma_',1,1,T);
df_F_1 = ones(T,1)*10;
df_F_2 = ones(T,1)*15;

options = optimoptions('fmincon','Display','final','MaxFunEval',1e5);
eparam = NaN(12,7);
parfor ii = 1:12
    data_mat = mncsFrnd( df_F_1, df_F_2, Sigma_, Omega_);

    fun_ = @(x) fun(x,data_mat);

    eparam(ii,:) = fmincon(fun_,[25,25,1,0,1,1,1],[],[],[],[],[3;3],[],[],options);
    eparam(ii,:) = fmincon(fun_,eparam(ii,:),[],[],[],[],[3;3],[],[],options);
end

eparam
for ii = 1:12
    eSigma(:,:,ii) = ivech(eparam(ii,3:5),'lower')*ivech(eparam(ii,3:5),'lower')';
    eOmega(:,:,ii) = eparam(ii,6:7)'*eparam(ii,6:7);
    eEV(:,:,ii) = eSigma(:,:,ii) + eOmega(:,:,ii)/eparam(ii,1);
end
%%
function y = fun(x,data_mat)
if size(x,1) < size(x,2)
    x = x';
end    
if size(x,1) ~= 7
    error('Size')
end
T = size(data_mat,3);
df_F_1 = ones(T,1)*x(1);
df_F_2 = ones(T,1)*x(2);
Sigma_ = ivech(x(3:5),'lower');
Sigma_ = repmat(Sigma_*Sigma_',1,1,T);
Omega_ = repmat(x(6:7)*x(6:7)',1,1,T);
y = mncsFlike(data_mat, df_F_1, df_F_2, Sigma_, Omega_, 100);
end
