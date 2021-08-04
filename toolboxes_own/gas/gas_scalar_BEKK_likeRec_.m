function [ nLogL, logLcontr, SigmaE, ScaledScore, nu, nuScore, varargout ] = ...
    gas_scalar_BEKK_likeRec_( param, p, q, r, s, X, dist )
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
elseif strcmp( dist, 'iWish' )
    n = param(k_+ p + q + 1);    
elseif strcmp( dist, 'tWish' )
    n = param(k_+ p + q + 1); 
    nu_intrcpt = param(k_+ p + q + 2);
    nu_scoreparam = param(k_+ p + q + 2 + r);
    nu_garchparam = param(k_+ p + q + 2 + r + s);
elseif strcmp( dist, 'itWish' )
    n = param(k_+ p + q + 1); 
    nu_intrcpt = param(k_+ p + q + 2);
    nu_scoreparam = param(k_+ p + q + 2 + r);
    nu_garchparam = param(k_+ p + q + 2 + r + s);    
elseif strcmp( dist, 'F' )
    n = param(k_+ p + q + 1); 
    nu_intrcpt = param(k_+ p + q + 2);
    nu_scoreparam = param(k_+ p + q + 2 + r);
    nu_garchparam = param(k_+ p + q + 2 + r + s);    
elseif strcmp( dist, 'Riesz' )
    n = param(k_+ p + q + 1 : k_+ p + q + k);    
elseif strcmp( dist, 'iRiesz' )
    n = param(k_+ p + q + 1 : k_+ p + q + k);  
elseif strcmp( dist, 'tRiesz' )
    n = param(k_+ p + q + 1 : k_+ p + q + k); 
    nu_intrcpt = param(k_+ p + q + k + 1);
    nu_scoreparam = param(k_+ p + q + k + 1 + r);
    nu_garchparam = param(k_+ p + q + k + 1 + r + s);    
elseif strcmp( dist, 'itRiesz' )
    n = param(k_+ p + q + 1 : k_+ p + q + k); 
    nu_intrcpt = param(k_+ p + q + k + 1);
    nu_scoreparam = param(k_+ p + q + k + 1 + r);
    nu_garchparam = param(k_+ p + q + k + 1 + r + s);       
elseif strcmp( dist, 'FRiesz' )
    n = param(k_+ p + q + 1 : k_+ p + q + k); 
    nu_intrcpt = param(k_+ p + q + k + 1 : k_+ p + q + k + k);    
    nu_scoreparam = param(k_+ p + q + k + k + r);
    nu_garchparam = param(k_+ p + q + k + k + r + s);       
end

if nargout >= 7
    param_out.intrcpt = intrcpt;
    param_out.scoreparam = scoreparam;
    param_out.garch_param = garchparam;
    param_out.n = n;
    
    if exist('nu_intrcpt','var')
        param_out.nu_intrcpt = nu_intrcpt;
        param_out.nu_scoreparam = nu_scoreparam;
        param_out.nu_garchparam = nu_garchparam;
    end

    param_out.all = param;
    
    varargout{1} = param_out;
end
%% Data Storage
SigmaE = NaN(k,k,T+22);
ScaledScore = NaN(k,k,T);
nu = NaN(T+22,1);
nuScore = NaN(T,1);
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
    SigmaE(:,:,tt) = intrcpt;
    for jj = 1:p
        if (tt-jj) <= 0
            SigmaE(:,:,tt) = SigmaE(:,:,tt) + scoreparam(jj)*zeros(k);
        else
            SigmaE(:,:,tt) = SigmaE(:,:,tt) + scoreparam(jj)*ScaledScore(:,:,tt-jj);
        end
    end
    for jj = 1:q
        if (tt-jj) <= 0
            SigmaE(:,:,tt) = SigmaE(:,:,tt) + garchparam(jj)*ini;
        else
            SigmaE(:,:,tt) = SigmaE(:,:,tt) + garchparam(jj)*SigmaE(:,:,tt-jj);
        end
    end
    
    nu(tt) = nu_intrcpt;
    for jj = 1:r
        if (tt-jj) <= 0
            nu(tt) = nu(tt) + nu_scoreparam(jj)*0;
        else
            nu(tt) = nu(tt) + nu_scoreparam(jj)*nuScore(tt-jj);
        end
    end
    for jj = 1:s
        if (tt-jj) <= 0
            nu(tt) = nu(tt) + nu_garchparam(jj)*nu_intrcpt;
        else
            nu(tt) = nu(tt) + nu_garchparam(jj)*nu(tt-jj);
        end
    end    

    % Likelihood Evaluation     
    try
   
        if exist('nu','var')
            [logLcontr(tt), score] = logpdf( dist, X(:,:,tt), SigmaE(:,:,tt), n, nu(tt) );
        else
            [logLcontr(tt), score] = logpdf( dist, X(:,:,tt), SigmaE(:,:,tt), n );
        end
        nuScore(tt) = score.nu;
        
        S = ivech(score.Sigma_);
        S = (S + diag(diag(S)))./2; %disregard symmetry in differential
    catch ME
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
    if strcmp( dist, 'Wish' )
        S = 2*SigmaE(:,:,tt)*S*SigmaE(:,:,tt);      
    elseif strcmp( dist, 'iWish' )
        S = -2/(n-k-1)/n*SigmaE(:,:,tt)*S*SigmaE(:,:,tt);
    elseif strcmp( dist, 'tWish')
        S = 2/n/nu(tt)*(nu(tt)-2)*SigmaE(:,:,tt)*S*SigmaE(:,:,tt);
    elseif strcmp( dist, 'itWish')
        S = -2/(n-k-1)/n/nu_intrcpt*(nu_intrcpt-2)*SigmaE(:,:,tt)*S*SigmaE(:,:,tt);
    elseif strcmp( dist, 'F' )
        S = 2/(n+1)*SigmaE(:,:,tt)*S*SigmaE(:,:,tt);
    elseif strcmp( dist, 'Riesz' )
        S = 2*SigmaE(:,:,tt)*S*SigmaE(:,:,tt);
    elseif strcmp( dist, 'iRiesz' )
        S = -2*k/sum(diag(M).^(-1))/k*sum(n.^(-1))*SigmaE(:,:,tt)*S*SigmaE(:,:,tt);
    elseif strcmp( dist, 'tRiesz' )
        S = 2*SigmaE(:,:,tt)*S*SigmaE(:,:,tt);
    elseif strcmp( dist, 'itRiesz' )
        S = -2*k/sum(diag(M).^(-1))*sum(n.^(-1))/k/nu_intrcpt*(nu_intrcpt-2)*SigmaE(:,:,tt)*S*SigmaE(:,:,tt);
    elseif strcmp( dist, 'FRiesz' )
        S = 2*SigmaE(:,:,tt)*S*SigmaE(:,:,tt);
    end
    ScaledScore(:,:,tt) = (S+S')./2; % to get rid of small asymmetries.
end
%% Fcst
for tt=T+1:T+22
    SigmaE(:,:,tt) = intrcpt;
    for jj = 1:p
        if (tt-jj) > T
            SigmaE(:,:,tt) = SigmaE(:,:,tt) + scoreparam(jj)*zeros(k);
        else
            SigmaE(:,:,tt) = SigmaE(:,:,tt) + scoreparam(jj)*ScaledScore(:,:,tt-jj);
        end
    end
    for jj = 1:q
        SigmaE(:,:,tt) = SigmaE(:,:,tt) + garchparam(jj)*SigmaE(:,:,tt-jj);
    end
    
    nu(tt) = nu_intrcpt;
    for jj = 1:r
        if (tt-jj) > T
            nu(tt) = nu(tt) + nu_scoreparam(jj)*0;
        else
            nu(tt) = nu(tt) + nu_scoreparam(jj)*nuScore(tt-jj);
        end
    end
    for jj = 1:s
        nu(tt) = nu(tt) + nu_garchparam(jj)*nu(tt-jj);
    end     
end
%% Log-Likelihood
nLogL = -sum(logLcontr);
%% Fit-Plot(s)
if nargout >= 8
    fitplot = figure("Visible",false,"WindowState",'fullscreen');
    y1 = sum(diag3d(X),2);
    y2 = sum(diag3d(SigmaE(:,:,1:T)),2);
    y3 = nu(1:T);
    plot(y1);
    hold on
    plot(y2,'LineWidth',2)
    text(T/2,max(y1)*3/4,strcat(num2str(sum(logLcontr))," | Score-Parameters:", num2str(scoreparam), " | Garch-Parameters:", num2str(garchparam), " | Matrix Euclidean Distance:", num2str(mean(matrix_euclidean_loss(X,SigmaE(:,:,1:T))))));    
    yyaxis right
    plot(y3)
    varargout{2} = fitplot;
end
