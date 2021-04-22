function [nLogL,logLcontr,covMat,fcst_rho_] = dcc_like(paramDCC,innov_a,condVar_a)
[T,N] = size(innov_a);
innov_a_std = innov_a./sqrt(condVar_a);
% DCC Correlation Likelihood
[rho_,fcst_rho_] = dcc_rec(paramDCC,innov_a_std);
covMat = NaN(N,N,T);
for t = 1:T
    covMat(:,:,t) = rho_(:,:,t).*(sqrt(condVar_a(t,:))'*sqrt(condVar_a(t,:)));
end
[nLogL,logLcontr]=mvnormlike(innov_a,zeros(size(innov_a)),covMat);
end