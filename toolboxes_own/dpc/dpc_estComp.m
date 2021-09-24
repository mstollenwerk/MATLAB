function [eparamComp] = dpc_estComp( covMdata, paramLoad, lP, lQ, cP, cQ, x0)
[n,~,T]=size(covMdata);
%% Get diagLRL
eQ = sbekkRec_s(paramLoad,covMdata,lP,lQ);
[eL_,~]=dpceig(mean(covMdata,3)); % Get Intercept Eigenvectors to enforce sign convention below. (Assumptions 2.1.2. - 2.1.3.)
diagLRL = NaN(n,T); % Get component input data
for t=1:T
    [eL,~] = dpceig(eQ{t},eL_); 
    diagLRL(:,t) = diag(eL'*covMdata(:,:,t)*eL);
end
%%
[ eparamComp, e_d, fcstd ] = dpc_estimatorComp( covMdata, diagLRL, [], cP, cQ, x0 );

end