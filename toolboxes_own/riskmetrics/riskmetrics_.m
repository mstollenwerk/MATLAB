function [Sigma_, fcst_Sigma] = riskmetrics_( data_RC )
%RISKMETRICS

[N,~,T] = size(data_RC);

Sigma_ = NaN(N,N,T+1);
Sigma_(:,:,1) = mean(data_RC,3);
for tt = 2:T+1
    Sigma_(:,:,tt) = .06*data_RC(:,:,tt-1) + .94*Sigma_(:,:,tt-1);
end

fcst_Sigma = Sigma_(:,:,T+1);
Sigma_ = Sigma_(:,:,1:T);

end