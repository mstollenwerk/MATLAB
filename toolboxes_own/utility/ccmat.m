function [ ROH_hat ] = ccmat( rt, lags )
%sCCM is Sample Cross-Correlation Matrix at various lags of Time Series R, 
%where rows of R are the different Time Series and columns of R are 
%different points in time

R=(rt-diag(mean(rt,2))*ones(size(rt)));
Rl=(rt-diag(mean(rt,2))*ones(size(rt)))';

sCCoMzero=(R*Rl)/size(R,2);

D_hat=sqrt(diag(diag(sCCoMzero)));

GAMMA_hat=(R(:,1:(size(R,2)-lags))*Rl((lags+1):size(Rl,1),:))/size(R,2);

ROH_hat=inv(D_hat)*GAMMA_hat*inv(D_hat);


end