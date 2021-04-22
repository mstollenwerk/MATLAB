function [rho_, fcst_rho_] = dcc_rec(paramDCC, innov_a_std)
% DCC Correlation Recursion
    [T,n]=size(innov_a_std);
    corr_mat = NaN(n,n,T);
    for t=1:T
        corr_mat(:,:,t) = innov_a_std(t,:)'*innov_a_std(t,:);
    end
    [ Q, fcst_Q ] = sbekkRec_s( paramDCC, corr_mat, 1, 1 );
    rho_ = NaN(n,n,T);
    for t=1:T
        [~,rho_(:,:,t)] = cov2corr(Q{t});
    end
    [~,fcst_rho_] = cov2corr(fcst_Q);
end