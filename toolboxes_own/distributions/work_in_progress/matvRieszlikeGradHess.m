function [ nLogL, g, H ] = matvRieszlikeGradHess( Omega_, n, R )

[p,~,N] = size(R);
%% Param
COm = chol(Omega_,'lower');
%% Log-likelihood computation
logLcontr = NaN(N,1);

term1 = -sum(n)/2*log(2);
term2 = -lgmvgammaln(n./2);
term3 = -loglpwdet([],n./2, diag(COm));

log_normalizing_constant = term1 + term2 + term3;

for ii = 1:N
    term4 = loglpwdet(R(:,:,ii),(n-p-1)./2);
    term5 = -trace(Omega_\R(:,:,ii))./2;
    
    log_kernel = term4 + term5;

    logLcontr(ii) = log_normalizing_constant + log_kernel;
end
nLogL = -sum(logLcontr);
%% Gradient and Hessian
    
gradient_ii = NaN(N,p);

for ii = 1:N

    gradient_ii(ii,:) = .5*( -log(2) - mvpsi(n/2) ) ...
                    + log(diag(chol(R(:,:,ii),'lower'))) - log(diag(COm));

end

g = -sum(gradient_ii,1)';

H = diag(N*.25*mvpsi(n/2,1));

end
