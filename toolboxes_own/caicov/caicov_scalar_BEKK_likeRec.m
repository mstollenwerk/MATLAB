function [ nLogL, logLcontr, Sigma_, varargout ] = ...
    caicov_scalar_BEKK_likeRec( param, p, q, X, dist )
%
% Michael Stollenwerk
% michael.stollenwerk@live.com
% 05.04.2021

[k,~,T] = size(X);
k_ = k*(k+1)/2;
%% Parameters
% Initialize with Backcast. Code copied from Sheppard MFE Toolbox.
m = ceil(sqrt(T));
w = .06 * .94.^(0:(m-1));
w = reshape(w/sum(w),[1 1 m]);
backCast = sum(bsxfun(@times,w,X(:,:,1:m)),3);
% % Initialize recursion at unconditional mean (stationarity assumed).
% ini = intrcpt./(1 - sum(garchparam) - sum(archparam));

intrcpt = ivechchol(param(1:k_));
archparam = param(k_ + 1 : k_ + p);
garchparam = param(k_ + p + 1 : k_+ p + q);

if strcmp( dist, 'Wish' )
    n = param(k_+ p + q + 1);
    ini = backCast/n;
    loglike = @(x1,x2,x3) matvWishlike(x1,x2,x3);
elseif strcmp( dist, 'iWish' )
    n = param(k_+ p + q + 1);
    ini = backCast*(n-k-1); 
    loglike = @(x1,x2,x3) matviWishlike(x1,x2,x3);   
elseif strcmp( dist, 'tWish' )
    n = param(k_+ p + q + 1); 
    nu = param(k_+ p + q + 2);
    ini = backCast/n/nu*(nu-2);  
    loglike = @(x1,x2,x3,x4) matvtWishlike(x1,x2,x3,x4);    
elseif strcmp( dist, 'itWish' )
    n = param(k_+ p + q + 1); 
    nu = param(k_+ p + q + 2);
    ini = backCast*(n-k-1); 
    loglike = @(x1,x2,x3,x4) matvitWishlike(x1,x2,x3,x4);
elseif strcmp( dist, 'F' )
    n = param(k_+ p + q + 1); 
    nu = param(k_+ p + q + 2);
    ini = backCast/n/nu*(n-k-1);
    loglike = @(x1,x2,x3,x4) matvFlike(x1,x2,x3,x4);    
elseif strcmp( dist, 'Riesz' )
    n = param(k_+ p + q + 1 : k_+ p + q + k);
    ini = zeros(k);
    for ii = 1:m
        CCC = chol(X(:,:,ii),'lower')/sqrtm(diag(n));
        ini = ini + w(ii)*(CCC*CCC');
    end
    loglike = @(x1,x2,x3) matvRieszlike(x1,x2,x3);    
elseif strcmp( dist, 'Riesz2' )
    n = param(k_+ p + q + 1 : k_+ p + q + k);
    for ii = 1:m
        CCC = cholU(X(:,:,ii),'lower')/sqrtm(diag(n));
        ini = ini + w(ii)*(CCC*CCC');
    end
    loglike = @(x1,x2,x3) matvRiesz2like(x1,x2,x3);        
elseif strcmp( dist, 'iRiesz' )
    n = param(k_+ p + q + 1 : k_+ p + q + k);  
    ini = zeros(k);
    for ii = 1:m
        CCC = chol(X(:,:,ii),'lower')/sqrtm(matviRieszexpmat(n));
        ini = ini + w(ii)*(CCC*CCC');
    end
    loglike = @(x1,x2,x3) matviRieszlike(x1,x2,x3); 
elseif strcmp( dist, 'iRiesz2' )
    n = param(k_+ p + q + 1 : k_+ p + q + k); 
    ini = zeros(k);
    for ii = 1:m
        CCC = chol(X(:,:,ii),'lower')/sqrtm(matviRiesz2expmat(n));
        ini = ini + w(ii)*(CCC*CCC');
    end
    loglike = @(x1,x2,x3) matviRiesz2like(x1,x2,x3); 
elseif strcmp( dist, 'tRiesz' )
    n = param(k_+ p + q + 1 : k_+ p + q + k); 
    nu = param(k_+ p + q + k + 1);
    ini = zeros(k);
    for ii = 1:m
        CCC = chol(X(:,:,ii),'lower')/sqrtm(diag(n)*nu/(nu-2));
        ini = ini + w(ii)*(CCC*CCC');
    end
    loglike = @(x1,x2,x3,x4) matvtRieszlike(x1,x2,x3,x4); 
elseif strcmp( dist, 'tRiesz2' )
    n = param(k_+ p + q + 1 : k_+ p + q + k); 
    nu = param(k_+ p + q + k + 1);
    ini = zeros(k);
    for ii = 1:m
        CCC = cholU(X(:,:,ii),'lower')/sqrtm(diag(n)*nu/(nu-2));
        ini = ini + w(ii)*(CCC*CCC');
    end    
    loglike = @(x1,x2,x3,x4) matvtRiesz2like(x1,x2,x3,x4); 
elseif strcmp( dist, 'itRiesz' )
    n = param(k_+ p + q + 1 : k_+ p + q + k); 
    nu = param(k_+ p + q + k + 1);
    ini = zeros(k);
    for ii = 1:m
        CCC = cholU(X(:,:,ii),'lower')/sqrtm(matviRieszexpmat(n));
        ini = ini + w(ii)*(CCC*CCC');
    end   
    loglike = @(x1,x2,x3,x4) matvitRieszlike(x1,x2,x3,x4); 
elseif strcmp( dist, 'itRiesz2' )
    n = param(k_+ p + q + 1 : k_+ p + q + k); 
    nu = param(k_+ p + q + k + 1);
    ini = zeros(k);
    for ii = 1:m
        CCC = chol(X(:,:,ii),'lower')/sqrtm(matviRiesz2expmat(n));
        ini = ini + w(ii)*(CCC*CCC');
    end
    loglike = @(x1,x2,x3,x4) matvitRiesz2like(x1,x2,x3,x4); 
elseif strcmp( dist, 'FRiesz' )
    n = param(k_+ p + q + 1 : k_+ p + q + k); 
    nu = param(k_+ p + q + k + 1 : k_+ p + q + k + k); 
    for ii = 1:m
        CCC = chol(X(:,:,ii),'lower')/sqrtm(matvFRieszexpmat(n,nu));
        ini = ini + w(ii)*(CCC*CCC');
    end
    loglike = @(x1,x2,x3,x4) matvFRieszlike(x1,x2,x3,x4); 
elseif strcmp( dist, 'FRiesz2' )
    n = param(k_+ p + q + 1 : k_+ p + q + k); 
    nu = param(k_+ p + q + k + 1 : k_+ p + q + k + k);  
    for ii = 1:m
        CCC = cholU(X(:,:,ii),'lower')/sqrtm(matvFRiesz2expmat(n,nu));
        ini = ini + w(ii)*(CCC*CCC');
    end
    loglike = @(x1,x2,x3,x4) matvFRiesz2like(x1,x2,x3,x4); 
end

if nargout >= 4
    param_out.intrcpt = intrcpt;
    param_out.archparam = archparam;
    param_out.garchparam = garchparam;
    param_out.n = n;
    
    if exist('nu','var')
        param_out.nu = nu;
    end

    param_out.all = param;
    
    varargout{1} = param_out;
end
%% Data Storage
Sigma_ = NaN(k,k,T+22);
logLcontr = NaN(T,1);
%% Recursion

for tt=1:T
    
    Sigma_(:,:,tt) = intrcpt;
    for jj = 1:p
        if (tt-jj) <= 0
            Sigma_(:,:,tt) = Sigma_(:,:,tt) + archparam(jj)*ini;
        else
            Sigma_(:,:,tt) = Sigma_(:,:,tt) + archparam(jj)*X(:,:,tt-jj);
        end
    end
    for jj = 1:q
        if (tt-jj) <= 0
            Sigma_(:,:,tt) = Sigma_(:,:,tt) + garchparam(jj)*ini;
        else
            Sigma_(:,:,tt) = Sigma_(:,:,tt) + garchparam(jj)*Sigma_(:,:,tt-jj);
        end
    end

    try
        % Likelihood Evaluation    
        if exist('nu','var')
            [nLogLcontr, ~, score, ~, ~, fisherinfo] = loglike( Sigma_(:,:,tt), n, nu, X(:,:,tt) );
        else
            [nLogLcontr, ~, score, ~, ~, fisherinfo] = loglike( Sigma_(:,:,tt), n, X(:,:,tt) );
        end
    logLcontr(tt) = -nLogLcontr;
    catch ME
%         tt
%         ME.message
%         [ME.stack.line]'
%         {ME.stack.name}'
        nLogL = inf;
        logLcontr = NaN;
        varargout{1} = NaN;      
        varargout{2} = NaN;
        return
    end
    
end
%% Fcst
for tt=T+1:T+22
    Sigma_(:,:,tt) = intrcpt;
    for jj = 1:p
        if (tt-jj) > T
            Sigma_(:,:,tt) = Sigma_(:,:,tt) + archparam(jj)*Sigma_(:,:,tt-jj);
        else
            Sigma_(:,:,tt) = Sigma_(:,:,tt) + archparam(jj)*X(:,:,tt-jj);
        end
    end
    for jj = 1:q
        Sigma_(:,:,tt) = Sigma_(:,:,tt) + garchparam(jj)*Sigma_(:,:,tt-jj);
    end
end
%% Log-Likelihood
nLogL = -sum(logLcontr);
%% Fit-Plot(s)
if nargout >= 5
    fitplot = NaN;
%     fitplot = figure("Visible",false,"WindowState",'fullscreen');
%     y1 = sum(diag3d(X),2);
%     y2 = sum(diag3d(SigmaE(:,:,1:T)),2);
%     plot(y1);
%     hold on
%     plot(y2,'LineWidth',2)
%     text(T/2,max(y1)*3/4,strcat(num2str(sum(logLcontr))," | ARCH-Parameters:", num2str(archparam), " | GARCH-Parameters:", num2str(garchparam), " | Matrix Euclidean Distance:", num2str(mean(matrix_euclidean_loss(X,SigmaE(:,:,1:T))))));
    varargout{2} = fitplot;
end
