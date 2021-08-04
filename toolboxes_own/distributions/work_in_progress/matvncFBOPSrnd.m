function Y = matvncFBOPSrnd(df_F_1, df_F_2, Sigma_, Omega_, N)
    
p = size(Sigma_,1);
Y = NaN(p,p,N);
for ii=1:N
    
    sqrtSig = sqrtm(Sigma_);
    V = matvncWishrnd(df_F_1, eye(p), (sqrtSig\Omega_)/sqrtSig, 1);
    T = wishrnd(eye(p),df_F_2);
    Y(:,:,ii) = (df_F_2 - p - 1)/df_F_1* sqrtSig*sqrtm(V)*inv(T)*sqrtm(V)*sqrtSig;
    
end

end