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
backCast = sum(bsxfun(@times,w,X(:,:,1:m)),3);
% % Initialize recursion at unconditional mean (stationarity assumed).
% ini = intrcpt./(1 - sum(garchparam) - sum(archparam));

intrcpt = ivechchol(param(1:k_));
archparam = param(k_ + 1 : k_ + p);
garchparam = param(k_ + p + 1 : k_+ p + q);

if nargout >= 4
    param_out.intrcpt = intrcpt;
    param_out.archparam = archparam;
    param_out.garchparam = garchparam;

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
            SigmaE(:,:,tt) = SigmaE(:,:,tt) + archparam(jj)*backCast;
        else
            SigmaE(:,:,tt) = SigmaE(:,:,tt) + archparam(jj)*X(:,:,tt-jj);
        end
    end
    for jj = 1:q
        if (tt-jj) <= 0
            SigmaE(:,:,tt) = SigmaE(:,:,tt) + garchparam(jj)*backCast;
        else
            SigmaE(:,:,tt) = SigmaE(:,:,tt) + garchparam(jj)*SigmaE(:,:,tt-jj);
        end
    end

    try
        % Centering distributions
        Sigma_ = matvStandardize(dist, SigmaE(:,:,tt), param(k_+p+q+1:end));
        % Likelihood Evaluation    
        nLogLcontr = matvLogLike( dist, Sigma_, param(k_+p+q+1:end), X(:,:,tt) );
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
