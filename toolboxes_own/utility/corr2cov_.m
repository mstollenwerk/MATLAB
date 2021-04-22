function covmat = corr2cov_(vars, corrs)
%% Input checking
[N,K] = size(vars);
[K_,K__,N_] = size(corrs);
if K~=K_ || N~=N_ || K_~=K__
    error('Input array sizes do not match')
end

%%
for ii = 1:N
    covmat(:,:,ii) = corrs(:,:,ii).*(sqrt(vars(ii,:)')*sqrt(vars(ii,:)));
end

end