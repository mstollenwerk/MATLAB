function [ nLogL, logLcontr, SigmaE, ScaledScore, varargout ] = ...
    gas_scalar_BEKK_likeRec( param, p, q, X, dist, scalingtype )
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

% Initialize with Backcast. Inspired by Sheppard MFE Toolbox.
m = ceil(sqrt(T));
w = .06 * .94.^(0:(m-1));
w = reshape(w/sum(w),[1 1 m]);
backCast_data = matvStandardize(dist, X(:,:,1:m), param(k_+p+q+1:end));
backCast = sum(bsxfun(@times,w,backCast_data(:,:,1:m)),3);
% Initialize recursion at unconditional mean (stationarity assumed).
% ini = intrcpt./(1 - sum(garchparam));
%% Data Storage
Sigma_ = NaN(k,k,T+22);
ScaledScore = NaN(k,k,T);
logLcontr = NaN(T,1);
%% Recursion
for tt=1:T
    
    Sigma_(:,:,tt) = intrcpt;
    for jj = 1:p
        if (tt-jj) <= 0
            Sigma_(:,:,tt) = Sigma_(:,:,tt) + scoreparam(jj)*zeros(k);
        else
            Sigma_(:,:,tt) = Sigma_(:,:,tt) + scoreparam(jj)*ScaledScore(:,:,tt-jj);
        end
    end
    for jj = 1:q
        if (tt-jj) <= 0
            Sigma_(:,:,tt) = Sigma_(:,:,tt) + garchparam(jj)*backCast;
        else
            Sigma_(:,:,tt) = Sigma_(:,:,tt) + garchparam(jj)*Sigma_(:,:,tt-jj);
        end
    end
    
    try
        % Likelihood Evaluation         
        [~, logLcontr(tt), score, ~, param_dist, fisherinfo] = ...
            matvLogLike( dist, Sigma_(:,:,tt), param(k_+p+q+1:end), X(:,:,tt) );   
    catch ME
%         tt
%         ME.message
%         [ME.stack.line]'
%         {ME.stack.name}'
        nLogL = inf;
        logLcontr = NaN;
        SigmaE = NaN;
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
        if (tt-jj) > T
            Sigma_(:,:,tt) = Sigma_(:,:,tt) + scoreparam(jj)*zeros(k);
        else
            Sigma_(:,:,tt) = Sigma_(:,:,tt) + scoreparam(jj)*ScaledScore(:,:,tt-jj);
        end
    end
    for jj = 1:q
        Sigma_(:,:,tt) = Sigma_(:,:,tt) + garchparam(jj)*Sigma_(:,:,tt-jj);
    end
end
SigmaE = matvEV(dist, Sigma_, param(k_+p+q+1:end));
%% Log-Likelihood
nLogL = -sum(logLcontr);
%% Parameter Output
if nargout >= 5
    param_out.intrcpt = intrcpt;
    param_out.scoreparam = scoreparam;
    param_out.garch_param = garchparam;
    param_out.n = param_dist.n;
    
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
