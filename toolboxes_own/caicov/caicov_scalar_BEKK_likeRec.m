function [ nLogL, logLcontr, SigmaE, varargout ] = ...
    caicov_scalar_BEKK_likeRec( param, p, q, X, dist )
%
% Michael Stollenwerk
% michael.stollenwerk@live.com
% 05.04.2021

t_ahead = 220;
[k,~,T] = size(X);
k_ = k*(k+1)/2;
%% Parameters
% Initialize with Backcast. Code copied from Sheppard MFE Toolbox.
m = ceil(sqrt(T));
w = .06 * .94.^(0:(m-1));
w = reshape(w/sum(w),[1 1 m]);
ini = sum(bsxfun(@times,w,X(:,:,1:m)),3);
% % Initialize recursion at unconditional mean (stationarity assumed).
% ini = intrcpt./(1 - sum(garchparam) - sum(archparam));

intrcpt = ivechchol(param(1:k_));
archparam = param(k_ + 1 : k_ + p);
garchparam = param(k_ + p + 1 : k_+ p + q);

if strcmp( dist, 'Wish' )
    n = param(k_+ p + q + 1);
    loglike = @(x1,x2,x3) matvWishlike(x1,x2,x3);
elseif strcmp( dist, 'iWish' )
    n = param(k_+ p + q + 1);
    loglike = @(x1,x2,x3) matviWishlike(x1,x2,x3);   
elseif strcmp( dist, 'tWish' )
    n = param(k_+ p + q + 1); 
    nu = param(k_+ p + q + 2);
    loglike = @(x1,x2,x3,x4) matvtWishlike(x1,x2,x3,x4);    
elseif strcmp( dist, 'itWish' )
    n = param(k_+ p + q + 1); 
    nu = param(k_+ p + q + 2);
    loglike = @(x1,x2,x3,x4) matvitWishlike(x1,x2,x3,x4);
elseif strcmp( dist, 'F' )
    n = param(k_+ p + q + 1); 
    nu = param(k_+ p + q + 2);
    loglike = @(x1,x2,x3,x4) matvFlike(x1,x2,x3,x4);    
elseif strcmp( dist, 'Riesz' )
    n = param(k_+ p + q + 1 : k_+ p + q + k);
    loglike = @(x1,x2,x3) matvRieszlike(x1,x2,x3);    
elseif strcmp( dist, 'Riesz2' )
    n = param(k_+ p + q + 1 : k_+ p + q + k);
    loglike = @(x1,x2,x3) matvRiesz2like(x1,x2,x3);        
elseif strcmp( dist, 'iRiesz' )
    n = param(k_+ p + q + 1 : k_+ p + q + k);  
    loglike = @(x1,x2,x3) matviRieszlike(x1,x2,x3); 
elseif strcmp( dist, 'iRiesz2' )
    n = param(k_+ p + q + 1 : k_+ p + q + k); 
    loglike = @(x1,x2,x3) matviRiesz2like(x1,x2,x3); 
elseif strcmp( dist, 'tRiesz' )
    n = param(k_+ p + q + 1 : k_+ p + q + k); 
    nu = param(k_+ p + q + k + 1);
    loglike = @(x1,x2,x3,x4) matvtRieszlike(x1,x2,x3,x4); 
elseif strcmp( dist, 'tRiesz2' )
    n = param(k_+ p + q + 1 : k_+ p + q + k); 
    nu = param(k_+ p + q + k + 1);
    loglike = @(x1,x2,x3,x4) matvtRiesz2like(x1,x2,x3,x4); 
elseif strcmp( dist, 'itRiesz' )
    n = param(k_+ p + q + 1 : k_+ p + q + k); 
    nu = param(k_+ p + q + k + 1);
    loglike = @(x1,x2,x3,x4) matvitRieszlike(x1,x2,x3,x4); 
elseif strcmp( dist, 'itRiesz2' )
    n = param(k_+ p + q + 1 : k_+ p + q + k); 
    nu = param(k_+ p + q + k + 1);
    loglike = @(x1,x2,x3,x4) matvitRiesz2like(x1,x2,x3,x4); 
elseif strcmp( dist, 'FRiesz' )
    n = param(k_+ p + q + 1 : k_+ p + q + k); 
    nu = param(k_+ p + q + k + 1 : k_+ p + q + k + k);
    loglike = @(x1,x2,x3,x4) matvFRieszlike(x1,x2,x3,x4); 
elseif strcmp( dist, 'FRiesz2' )
    n = param(k_+ p + q + 1 : k_+ p + q + k); 
    nu = param(k_+ p + q + k + 1 : k_+ p + q + k + k);
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
SigmaE = NaN(k,k,T+t_ahead);
logLcontr = NaN(T,1);
%% Recursion

for tt=1:T
    
    SigmaE(:,:,tt) = intrcpt;
    for jj = 1:p
        if (tt-jj) <= 0
            SigmaE(:,:,tt) = SigmaE(:,:,tt) + archparam(jj)*ini;
        else
            SigmaE(:,:,tt) = SigmaE(:,:,tt) + archparam(jj)*X(:,:,tt-jj);
        end
    end
    for jj = 1:q
        if (tt-jj) <= 0
            SigmaE(:,:,tt) = SigmaE(:,:,tt) + garchparam(jj)*ini;
        else
            SigmaE(:,:,tt) = SigmaE(:,:,tt) + garchparam(jj)*SigmaE(:,:,tt-jj);
        end
    end

    try
        % Centering distributions
        if strcmp( dist, 'Wish' )
            Sigma_ = SigmaE(:,:,tt)/n;
        elseif strcmp( dist, 'iWish' )
            Sigma_ = SigmaE(:,:,tt)*(n-k-1);
        elseif strcmp( dist, 'tWish')
            Sigma_ = SigmaE(:,:,tt)*(nu-2)/nu/n;
        elseif strcmp( dist, 'itWish')
            Sigma_ = SigmaE(:,:,tt)*(n-k-1);
        elseif strcmp( dist, 'F' )
            Sigma_ = SigmaE(:,:,tt)*(nu-k-1)/n/nu;
        elseif strcmp( dist, 'Riesz' )
            C = chol(SigmaE(:,:,tt),'lower');
            Sigma_ = C/diag(n)*C';
        elseif strcmp( dist, 'Riesz2' )
            C = cholU(SigmaE(:,:,tt));
            Sigma_ = C/diag(n)*C';            
        elseif strcmp( dist, 'iRiesz' )
            C = cholU(SigmaE(:,:,tt));
            M = matviRieszexpmat(n);
            Sigma_ = C/M*C';
        elseif strcmp( dist, 'iRiesz2' )
            C = chol(SigmaE(:,:,tt),'lower');
            Sigma_ = C/matviRiesz2expmat(n)*C';            
        elseif strcmp( dist, 'tRiesz' )
            C = chol(SigmaE(:,:,tt),'lower');
            M = matvtRieszexpmat(n,nu);
            Sigma_ = C/M*C';
        elseif strcmp( dist, 'tRiesz2' )
            C = cholU(SigmaE(:,:,tt));
            M = matvtRieszexpmat(n,nu);
            Sigma_ = C/M*C';            
        elseif strcmp( dist, 'itRiesz' )
            C = cholU(SigmaE(:,:,tt));
            M = matviRieszexpmat(n);
            Sigma_ = C/M*C';
        elseif strcmp( dist, 'itRiesz2' )
            C = chol(SigmaE(:,:,tt),'lower');
            M = matviRieszexpmat(n);
            Sigma_ = C/M*C';            
        elseif strcmp( dist, 'FRiesz' )
            C = chol(SigmaE(:,:,tt),'lower');
            M = matvFRieszexpmat(n,nu);
            Sigma_ = C/M*C';
        elseif strcmp( dist, 'FRiesz2' )
            C = cholU(SigmaE(:,:,tt));
            M = matvFRiesz2expmat(n,nu);
            Sigma_ = C/M*C';            
        end
        % Likelihood Evaluation    
        if exist('nu','var')
            nLogLcontr = loglike( Sigma_, n, nu, X(:,:,tt) );
        else
            nLogLcontr = loglike( Sigma_, n, X(:,:,tt) );
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
for tt=T+1:T+t_ahead
    SigmaE(:,:,tt) = intrcpt;
    for jj = 1:p
        if (tt-jj) > T
            SigmaE(:,:,tt) = SigmaE(:,:,tt) + archparam(jj)*SigmaE(:,:,tt-jj);
        else
            SigmaE(:,:,tt) = SigmaE(:,:,tt) + archparam(jj)*X(:,:,tt-jj);
        end
    end
    for jj = 1:q
        SigmaE(:,:,tt) = SigmaE(:,:,tt) + garchparam(jj)*SigmaE(:,:,tt-jj);
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
