function [ nLogL, logLcontr, Sigma_, ScaledScore, varargout ] = ...
    gas_diag_BEKK_likeRec( param, p, q, X, dist, scalingtype )
%
% Michael Stollenwerk
% michael.stollenwerk@live.com
% 24.09.2021

[k,~,T] = size(X);
k_ = k*(k+1)/2;
%% Parameters
% Initialize with Backcast. Code copied from Sheppard MFE Toolbox.
m = ceil(sqrt(T));
w = .06 * .94.^(0:(m-1));
w = reshape(w/sum(w),[1 1 m]);
backCast = sum(bsxfun(@times,w,X(:,:,1:m)),3);
% Initialize recursion at unconditional mean (stationarity assumed).
% ini = intrcpt./(1 - sum(garchparam));

intrcpt = ivechchol(param(1:k_));
scoreparam = reshape(param(k_ + 1 : k_ + p*k),k,p);
garchparam = reshape(param(k_ + p*k + 1 : k_+ p*k + q*k),k,q);

if strcmp( dist, 'Wish' )
    n = param(k_+ p*k + q*k + 1);
    ini = backCast/n;
    loglike = @(x1,x2,x3) matvWishlike(x1,x2,x3);
elseif strcmp( dist, 'iWish' )
    n = param(k_+ p*k + q*k + 1);
    ini = backCast*(n-k-1); 
    loglike = @(x1,x2,x3) matviWishlike(x1,x2,x3);   
elseif strcmp( dist, 'tWish' )
    n = param(k_+ p*k + q*k + 1); 
    nu = param(k_+ p*k + q*k + 2);
    ini = backCast/n/nu*(nu-2);  
    loglike = @(x1,x2,x3,x4) matvtWishlike(x1,x2,x3,x4);    
elseif strcmp( dist, 'itWish' )
    n = param(k_+ p*k + q*k + 1); 
    nu = param(k_+ p*k + q*k + 2);
    ini = backCast*(n-k-1); 
    loglike = @(x1,x2,x3,x4) matvitWishlike(x1,x2,x3,x4);
elseif strcmp( dist, 'F' )
    n = param(k_+ p*k + q*k + 1); 
    nu = param(k_+ p*k + q*k + 2);
    ini = backCast/n/nu*(nu-k-1);
    loglike = @(x1,x2,x3,x4) matvFlike(x1,x2,x3,x4);    
elseif strcmp( dist, 'Riesz' )
    n = param(k_+ p*k + q*k + 1 : k_+ p*k + q*k + k);
    ini = zeros(k);
    for ii = 1:m
        CCC = chol(X(:,:,ii),'lower')/sqrtm(diag(n));
        ini = ini + w(ii)*(CCC*CCC');
    end
    loglike = @(x1,x2,x3) matvRieszlike(x1,x2,x3);    
elseif strcmp( dist, 'Riesz2' )
    n = param(k_+ p*k + q*k + 1 : k_+ p*k + q*k + k);
    ini = zeros(k);
    for ii = 1:m
        CCC = cholU(X(:,:,ii))/sqrtm(diag(n));
        ini = ini + w(ii)*(CCC*CCC');
    end
    loglike = @(x1,x2,x3) matvRiesz2like(x1,x2,x3);        
elseif strcmp( dist, 'iRiesz' )
    n = param(k_+ p*k + q*k + 1 : k_+ p*k + q*k + k);  
    ini = zeros(k);
    for ii = 1:m
        CCC = chol(X(:,:,ii),'lower')/sqrtm(matviRieszexpmat(n));
        ini = ini + w(ii)*(CCC*CCC');
    end
    loglike = @(x1,x2,x3) matviRieszlike(x1,x2,x3); 
elseif strcmp( dist, 'iRiesz2' )
    n = param(k_+ p*k + q*k + 1 : k_+ p*k + q*k + k); 
    ini = zeros(k);
    for ii = 1:m
        CCC = chol(X(:,:,ii),'lower')/sqrtm(matviRiesz2expmat(n));
        ini = ini + w(ii)*(CCC*CCC');
    end
    loglike = @(x1,x2,x3) matviRiesz2like(x1,x2,x3); 
elseif strcmp( dist, 'tRiesz' )
    n = param(k_+ p*k + q*k + 1 : k_+ p*k + q*k + k); 
    nu = param(k_+ p*k + q*k + k + 1);
    ini = zeros(k);
    for ii = 1:m
        CCC = chol(X(:,:,ii),'lower')/sqrtm(diag(n)*nu/(nu-2));
        ini = ini + w(ii)*(CCC*CCC');
    end
    loglike = @(x1,x2,x3,x4) matvtRieszlike(x1,x2,x3,x4); 
elseif strcmp( dist, 'tRiesz2' )
    n = param(k_+ p*k + q*k + 1 : k_+ p*k + q*k + k); 
    nu = param(k_+ p*k + q*k + k + 1);
    ini = zeros(k);
    for ii = 1:m
        CCC = cholU(X(:,:,ii))/sqrtm(diag(n)*nu/(nu-2));
        ini = ini + w(ii)*(CCC*CCC');
    end    
    loglike = @(x1,x2,x3,x4) matvtRiesz2like(x1,x2,x3,x4); 
elseif strcmp( dist, 'itRiesz' )
    n = param(k_+ p*k + q*k + 1 : k_+ p*k + q*k + k); 
    nu = param(k_+ p*k + q*k + k + 1);
    ini = zeros(k);
    for ii = 1:m
        CCC = cholU(X(:,:,ii))/sqrtm(matviRieszexpmat(n));
        ini = ini + w(ii)*(CCC*CCC');
    end   
    loglike = @(x1,x2,x3,x4) matvitRieszlike(x1,x2,x3,x4); 
elseif strcmp( dist, 'itRiesz2' )
    n = param(k_+ p*k + q*k + 1 : k_+ p*k + q*k + k); 
    nu = param(k_+ p*k + q*k + k + 1);
    ini = zeros(k);
    for ii = 1:m
        CCC = chol(X(:,:,ii),'lower')/sqrtm(matviRiesz2expmat(n));
        ini = ini + w(ii)*(CCC*CCC');
    end
    loglike = @(x1,x2,x3,x4) matvitRiesz2like(x1,x2,x3,x4); 
elseif strcmp( dist, 'FRiesz' )
    n = param(k_+ p*k + q*k + 1 : k_+ p*k + q*k + k); 
    nu = param(k_+ p*k + q*k + k + 1 : k_+ p*k + q*k + k + k); 
    ini = zeros(k);    
    for ii = 1:m
        CCC = chol(X(:,:,ii),'lower')/sqrtm(matvFRieszexpmat(n,nu));
        ini = ini + w(ii)*(CCC*CCC');
    end
    loglike = @(x1,x2,x3,x4) matvFRieszlike(x1,x2,x3,x4); 
elseif strcmp( dist, 'FRiesz2' )
    n = param(k_+ p*k + q*k + 1 : k_+ p*k + q*k + k); 
    nu = param(k_+ p*k + q*k + k + 1 : k_+ p*k + q*k + k + k);  
    ini = zeros(k);    
    for ii = 1:m
        CCC = cholU(X(:,:,ii))/sqrtm(matvFRiesz2expmat(n,nu));
        ini = ini + w(ii)*(CCC*CCC');
    end
    loglike = @(x1,x2,x3,x4) matvFRiesz2like(x1,x2,x3,x4); 
end

if nargout >= 5
    param_out.intrcpt = intrcpt;
    param_out.scoreparam = scoreparam;
    param_out.garch_param = garchparam;
    param_out.n = n;
    
    if exist('nu','var')
        param_out.nu = nu;
    end

    param_out.all = param;
    
    varargout{1} = param_out;
end
%% Data Storage
Sigma_ = NaN(k,k,T+22);
ScaledScore = NaN(k,k,T);
logLcontr = NaN(T,1);
%% Recursion
for tt=1:T
    
    Sigma_(:,:,tt) = intrcpt;
    for jj = 1:p
        A = diag(sqrt(scoreparam(:,jj)));
        if (tt-jj) <= 0
            Sigma_(:,:,tt) = Sigma_(:,:,tt) + A*zeros(k)*A;
        else
            Sigma_(:,:,tt) = Sigma_(:,:,tt) + A*ScaledScore(:,:,tt-jj)*A;
        end
    end
    for jj = 1:q
        B = diag(sqrt(garchparam(:,jj)));
        if (tt-jj) <= 0
            Sigma_(:,:,tt) = Sigma_(:,:,tt) + B*ini*B;
        else
            Sigma_(:,:,tt) = Sigma_(:,:,tt) + B*Sigma_(:,:,tt-jj)*B;
        end
    end
    
    try
        % Likelihood Evaluation    
        if exist('nu','var')
            if scalingtype == 1
                [nLogLcontr, ~, score, ~, ~, fisherinfo] = loglike( Sigma_(:,:,tt), n, nu, X(:,:,tt) );
            elseif scalingtype == 2
                [nLogLcontr, ~, score ] = loglike( Sigma_(:,:,tt), n, nu, X(:,:,tt) );
            else
                error('Invalid scalingtype input. Must be 1 or 2.')
            end
        else
            if scalingtype == 1
                [nLogLcontr, ~, score, ~, ~, fisherinfo] = loglike( Sigma_(:,:,tt), n, X(:,:,tt) );
            elseif scalingtype == 2
                [nLogLcontr, ~, score ] = loglike( Sigma_(:,:,tt), n, X(:,:,tt) );
            else
                error('Invalid scalingtype input. Must be 1 or 2.')
            end
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
    
    % Scaled Score 
    if scalingtype == 1
        ScaledScore(:,:,tt) = ivech(fisherinfo.Sigma_\score.Sigma_');
    elseif scalingtype == 2
        ScaledScore(:,:,tt) = score.Sigma_WishFishScaling;
    else
        error('Invalid scalingtype input. Must be 1 or 2.')
    end
    
end
%% Fcst
for tt=T+1:T+22
    Sigma_(:,:,tt) = intrcpt;
    for jj = 1:p
        A = diag(sqrt(scoreparam(:,jj)));
        if (tt-jj) > T
            Sigma_(:,:,tt) = Sigma_(:,:,tt) + A*zeros(k)*A;
        else
            Sigma_(:,:,tt) = Sigma_(:,:,tt) + A*ScaledScore(:,:,tt-jj)*A;
        end
    end
    for jj = 1:q
        B = diag(sqrt(garchparam(:,jj)));
        Sigma_(:,:,tt) = Sigma_(:,:,tt) + B*Sigma_(:,:,tt-jj)*B;
    end
end
%% Log-Likelihood
nLogL = -sum(logLcontr);
