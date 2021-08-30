function [ nLogL, logLcontr, Omega_, ScaledScore, varargout ] = ...
    gas_scalar_BEKK_likeRec( param, p, q, X, dist )
%
% Michael Stollenwerk
% michael.stollenwerk@live.com
% 11.08.2020

[k,~,T] = size(X);
k_ = k*(k+1)/2;
%% Parameters
intrcpt = ivechchol(param(1:k_));
scoreparam = param(k_ + 1 : k_ + p);
garchparam = param(k_ + p + 1 : k_+ p + q);

if strcmp( dist, 'Wish' )
    n = param(k_+ p + q + 1);
    loglike = @(x1,x2,x3) matvsWishlike(x1,x2,x3);
elseif strcmp( dist, 'iWish' )
    n = param(k_+ p + q + 1);
    loglike = @(x1,x2,x3) matvsiWishlike(x1,x2,x3);   
elseif strcmp( dist, 'tWish' )
    n = param(k_+ p + q + 1); 
    nu = param(k_+ p + q + 2);
    loglike = @(x1,x2,x3,x4) matvstWishlike(x1,x2,x3,x4);    
elseif strcmp( dist, 'itWish' )
    n = param(k_+ p + q + 1); 
    nu = param(k_+ p + q + 2);
    loglike = @(x1,x2,x3,x4) matvsitWishlike(x1,x2,x3,x4);
elseif strcmp( dist, 'F' )
    n = param(k_+ p + q + 1); 
    nu = param(k_+ p + q + 2);
    loglike = @(x1,x2,x3,x4) matvsFlike(x1,x2,x3,x4);    
elseif strcmp( dist, 'Riesz' )
    n = param(k_+ p + q + 1 : k_+ p + q + k);
    loglike = @(x1,x2,x3) matvsRieszlike(x1,x2,x3);    
elseif strcmp( dist, 'iRiesz' )
    n = param(k_+ p + q + 1 : k_+ p + q + k);  
    loglike = @(x1,x2,x3) matvsiRieszlike(x1,x2,x3);    
elseif strcmp( dist, 'tRiesz' )
    n = param(k_+ p + q + 1 : k_+ p + q + k); 
    nu = param(k_+ p + q + k + 1);
    loglike = @(x1,x2,x3,x4) matvstRieszlike(x1,x2,x3,x4);    
elseif strcmp( dist, 'itRiesz' )
    n = param(k_+ p + q + 1 : k_+ p + q + k); 
    nu = param(k_+ p + q + k + 1);
    loglike = @(x1,x2,x3,x4) matvsitRieszlike(x1,x2,x3,x4); 
elseif strcmp( dist, 'itRiesz2' )
    n = param(k_+ p + q + 1 : k_+ p + q + k); 
    nu = param(k_+ p + q + k + 1);
    loglike = @(x1,x2,x3,x4) matvsitRiesz2like(x1,x2,x3,x4); 
elseif strcmp( dist, 'FRiesz' )
    n = param(k_+ p + q + 1 : k_+ p + q + k); 
    nu = param(k_+ p + q + k + 1 : k_+ p + q + k + k);        
    loglike = @(x1,x2,x3,x4) matvsFRieszlike(x1,x2,x3,x4);    
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
Omega_ = NaN(k,k,T+22);
ScaledScore = NaN(k,k,T);
logLcontr = NaN(T,1);
%% Recursion
% Initialize with Backcast. Code copied from Sheppard MFE Toolbox.
m = ceil(sqrt(T));
w = .06 * .94.^(0:(m-1));
w = reshape(w/sum(w),[1 1 m]);
backCast = sum(bsxfun(@times,w,X(:,:,1:m)),3);
ini = backCast;
% Initialize recursion at unconditional mean (stationarity assumed).
% ini = intrcpt./(1 - sum(garchparam));
for tt=1:T
    
    Omega_(:,:,tt) = intrcpt;
    for jj = 1:p
        if (tt-jj) <= 0
            Omega_(:,:,tt) = Omega_(:,:,tt) + scoreparam(jj)*zeros(k);
        else
            Omega_(:,:,tt) = Omega_(:,:,tt) + scoreparam(jj)*ScaledScore(:,:,tt-jj);
        end
    end
    for jj = 1:q
        if (tt-jj) <= 0
            Omega_(:,:,tt) = Omega_(:,:,tt) + garchparam(jj)*ini;
        else
            Omega_(:,:,tt) = Omega_(:,:,tt) + garchparam(jj)*Omega_(:,:,tt-jj);
        end
    end
    
    try
        % Likelihood Evaluation    
        if exist('nu','var')
            [logLcontr(tt), ~, score] = loglike( Omega_(:,:,tt), n, nu, X(:,:,tt) );
        else
            [logLcontr(tt), ~, score] = loglike( Omega_(:,:,tt), n, X(:,:,tt) );
        end
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
    ScaledScore(:,:,tt) = score.Omega_scaledbyiFish;
    
end
%% Fcst
for tt=T+1:T+22
    Omega_(:,:,tt) = intrcpt;
    for jj = 1:p
        if (tt-jj) > T
            Omega_(:,:,tt) = Omega_(:,:,tt) + scoreparam(jj)*zeros(k);
        else
            Omega_(:,:,tt) = Omega_(:,:,tt) + scoreparam(jj)*ScaledScore(:,:,tt-jj);
        end
    end
    for jj = 1:q
        Omega_(:,:,tt) = Omega_(:,:,tt) + garchparam(jj)*Omega_(:,:,tt-jj);
    end
end
%% Log-Likelihood
nLogL = -sum(logLcontr);
%% Fit-Plot(s)
if nargout >= 6
    fitplot = NaN;
%     fitplot = figure("Visible",false,"WindowState",'fullscreen');
%     y1 = sum(diag3d(X),2);
%     y2 = sum(diag3d(SigmaE(:,:,1:T)),2);
%     plot(y1);
%     hold on
%     plot(y2,'LineWidth',2)
%     text(T/2,max(y1)*3/4,strcat(num2str(sum(logLcontr))," | Score-Parameters:", num2str(scoreparam), " | Garch-Parameters:", num2str(garchparam), " | Matrix Euclidean Distance:", num2str(mean(matrix_euclidean_loss(X,SigmaE(:,:,1:T))))));
    varargout{2} = fitplot;
end
