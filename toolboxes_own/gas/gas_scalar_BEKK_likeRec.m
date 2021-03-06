function [ nLogL, logLcontr, Sigma_, S, varargout ] = ...
    gas_scalar_BEKK_likeRec( param, p, q, R, dist )
%
% Michael Stollenwerk
% michael.stollenwerk@live.com
% 11.08.2020

t_ahead = 220;

[k,~,T] = size(R);
k_ = k*(k+1)/2;
%% Parameters
intrcpt = ivechchol(param(1:k_));
scoreparam = reshape(param(k_ + 1 : k_ + 2*p),[],2);
garchparam = param(k_ + 2*p + 1 : k_+ 2*p + q);

% Initialize with Backcast. Inspired by Sheppard MFE Toolbox.
m = ceil(sqrt(T));
w = .06 * .94.^(0:(m-1));
w = reshape(w/sum(w),[1 1 m]);
backCast = sum(bsxfun(@times,w, R(:,:,1:m)),3);
iniSig = backCast;
% Initialize recursion at unconditional mean (stationarity assumed).
% iniSig = intrcpt./(1 - sum(garchparam));
%% Data Storage
Sigma_ = NaN(k,k,T+t_ahead);
S = NaN(k,k,T,2);
logLcontr = NaN(T,1);
%% Recursion
for tt=1:T
    
    Sigma_(:,:,tt) = intrcpt;
    for jj = 1:p
        if (tt-jj) <= 0
            Sigma_(:,:,tt) = Sigma_(:,:,tt);
        else
            Sigma_(:,:,tt) = Sigma_(:,:,tt) ...
                                + scoreparam(jj,1)*S(:,:,tt-jj,1) ...
                                + scoreparam(jj,2)*S(:,:,tt-jj,2);
        end
    end
    for jj = 1:q
        if (tt-jj) <= 0
            Sigma_(:,:,tt) = Sigma_(:,:,tt) + garchparam(jj)*iniSig;
        else
            Sigma_(:,:,tt) = Sigma_(:,:,tt) + garchparam(jj)*Sigma_(:,:,tt-jj);
        end
    end
    
    try
        % Likelihood Evaluation         
        [~, logLcontr(tt), score, ~, param_dist] = ...
            matvsLogLike( dist, Sigma_(:,:,tt), param(k_+2*p+q+1:end), R(:,:,tt) );   
    catch ME
%         tt
%         ME.message
%         [ME.stack.line]'
%         {ME.stack.name}'
        nLogL = inf;
        logLcontr = NaN;
        Sigma_ = NaN;
        varargout{1} = NaN;      
        varargout{2} = NaN;
        return
    end
    
    % Scaled Score 
    SigDel = Sigma_(:,:,tt)*score.SigmaNonSym;
    SigDelSig = SigDel*Sigma_(:,:,tt);
    S(:,:,tt,1) = SigDelSig+SigDelSig';
    S(:,:,tt,2) = trace(SigDel)*Sigma_(:,:,tt);
    
end
%% Fcst
for tt=T+1:T+t_ahead
    Sigma_(:,:,tt) = intrcpt;
    for jj = 1:p
        if (tt-jj) <= T

            Sigma_(:,:,tt) = Sigma_(:,:,tt) ...
                                + scoreparam(jj,1)*S(:,:,tt-jj,1) ...
                                + scoreparam(jj,2)*S(:,:,tt-jj,2);
        end
    end
    for jj = 1:q
        Sigma_(:,:,tt) = Sigma_(:,:,tt) + garchparam(jj)*Sigma_(:,:,tt-jj);
    end
% The commented code below is problematic, because it causes my_fmincon to 
% fail sometimes. If the true optimal likelihood lies at a point where the
% forecasts are non-pd, it tries to go closer and closer to this point,
% sometimes stopping at a point where there are non-finite derivatives,
% which causes my_fmincon to fail.
    if any(eig(Sigma_(:,:,tt))<0) %This is if the forecasted Sigmas are not pd, which can happen despite all in-sample sigmas being pd. in this case likely the estimated scoreparam is too high for stationarity.
        nLogL = inf;
        logLcontr = NaN;
        Sigma_ = NaN;
        varargout{1} = NaN;      
        varargout{2} = NaN;
        return
    end
end
% [Sigma_(:,:,T+1:T+t_ahead), pdadjustments_made] = makepd(Sigma_(:,:,T+1:T+t_ahead));
%% Log-Likelihood
nLogL = -sum(logLcontr);
%% Parameter Output
if nargout >= 5
    param_out.intrcpt = intrcpt;
    param_out.scoreparam = scoreparam;
    param_out.garch_param = garchparam;
    
    if isfield(param_dist,'n')
        param_out.n = param_dist.n;
    end
    
    if isfield(param_dist,'nu')
        param_out.nu = param_dist.nu;
    end

    param_out.all = param;
    
    varargout{1} = param_out;
end
%% Fit-Plot(s)
if nargout >= 6
    fitplot = NaN;
%     fitplot = figure("Visible",false,"WindowState",'fullscreen');
%     y1 = sum(diag3d(X),2);
%     for tt = 1:T
%         y2(tt) = sum(diag(cholU(Sigma_(:,:,tt))*matviRieszexpmat(n)*cholU(Sigma_(:,:,tt))'));
%     end
% %     y2 = sum(diag3d(SigmaE(:,:,1:T)),2);
%     plot(y1);
%     hold on
%     plot(y2,'LineWidth',2)
%     text(T/2,max(y1)*3/4,strcat(num2str(sum(logLcontr))," | Score-Parameters:", num2str(scoreparam), " | Garch-Parameters:", num2str(garchparam)));%, " | Matrix Euclidean Distance:", num2str(mean(matrix_euclidean_loss(X,SigmaE(:,:,1:T))))));
    varargout{2} = fitplot;
end
