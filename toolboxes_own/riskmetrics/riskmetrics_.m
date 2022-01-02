function [Sigma_, fcst_Sigma] = riskmetrics_( data_RC )
%RISKMETRICS

[N,~,T] = size(data_RC);

m = ceil(sqrt(T));
w = .06 * .94.^(0:(m-1));
w = reshape(w/sum(w),[1 1 m]);
backCast = sum(bsxfun(@times,w,data_RC(:,:,1:m)),3);

Sigma_ = NaN(N,N,T+1);
Sigma_(:,:,1) = backCast;
for tt = 2:T+1
    Sigma_(:,:,tt) = .06*data_RC(:,:,tt-1) + .94*Sigma_(:,:,tt-1);
end

fcst_Sigma = Sigma_(:,:,T+1);
Sigma_ = Sigma_(:,:,1:T);

end