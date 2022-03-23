function [ nLogL, logLcontr, Sigma_, ScaledScore, weights_out, varargout ] = ...
    gas_scalar_BEKK_mix_likeRec( param, p, q, R )
%
% Michael Stollenwerk
% michael.stollenwerk@live.com
% 11.08.2020

t_ahead = 220;

[k,~,T] = size(R);
k_ = k*(k+1)/2;
%% Parameters
intrcpt = ivechchol(param(1:k_));
scoreparam = param(k_ + 1 : k_ + p);
garchparam = param(k_ + p + 1 : k_+ p + q);
weightparam = param(k_ + p + q + 1 : k_ + p + q + 4);
n_itRiesz2 = param(k_ + p + q + 5 : k_ + p + q + 4 + k);
nu_itRiesz2 = param(k_ + p + q + 4 + k + 1);
n_FRiesz = param(k_ + p + q + 4 + k + 2 : k_ + p + q + 4 + k + 1 + k);
nu_FRiesz = param(k_ + p + q + 4 + k + 1 + k + 1: k_ + p + q + 4 + k + 1 + k + k);

% Initialize with Backcast. Inspired by Sheppard MFE Toolbox.
m = ceil(sqrt(T));
w = .06 * .94.^(0:(m-1));
w = reshape(w/sum(w),[1 1 m]);
backCast = sum(bsxfun(@times,w, R(:,:,1:m)),3);
% Initialize recursion at unconditional mean (stationarity assumed).
% ini = intrcpt./(1 - sum(garchparam));
%% Data Storage
Sigma_ = NaN(k,k,T+t_ahead);
ScaledScore = NaN(k,k,T);
weights = NaN(T+t_ahead,1);
logLcontr = NaN(T,1);
Lcontr_itRiesz2 = NaN(T,1);
Lcontr_FRiesz = NaN(T,1);
score_weights = NaN(T+t_ahead,1);
%% Recursion
weights(1) = weightparam(1);
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
    if tt>1
        weights(tt) = (1-weightparam(3))*weightparam(2) ...
                        + weightparam(3)*weights(tt-1) ...
                        + weightparam(4)*score_weights(tt-1);
    end   
    
    try
        % Likelihood Evaluation
        [nLogLcontr_itRiesz2, ~, scoreSig_itR] = matvsitRiesz2like(Sigma_(:,:,tt),n_itRiesz2,nu_itRiesz2,R(:,:,tt));
        [nLogLcontr_FRiesz, ~, scoreSig_FR] = matvsFRieszlike(Sigma_(:,:,tt),n_FRiesz,nu_FRiesz,R(:,:,tt));
        Lcontr_itRiesz2(tt) = exp(-nLogLcontr_itRiesz2);
        Lcontr_FRiesz(tt) = exp(-nLogLcontr_FRiesz);
        Like_tt = weights(tt)*Lcontr_itRiesz2(tt) + (1-weights(tt))*Lcontr_FRiesz(tt);
        score_weights(tt) = (Lcontr_itRiesz2(tt) - Lcontr_FRiesz(tt))/Like_tt;
        logLcontr(tt) = log(Like_tt);
    catch ME
%         tt
%         ME.message
%         [ME.stack.line]'
%         {ME.stack.name}'
        nLogL = inf;
        logLcontr = NaN;
        Sigma_ = NaN;
        ScaledScore = NaN;
        weights_out = NaN;
        varargout{1} = NaN;      
        varargout{2} = NaN;
        return
    end
    
    % Scaled Score 
    ScaledScore(:,:,tt) = weights(tt)*scoreSig_itR.rc_paper ...
                          + (1-weights(tt))*scoreSig_FR.rc_paper;
    
end
%% Fcst
for tt=T+1:T+t_ahead
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
    if tt == T+1
        weights(tt) = (1-weightparam(3))*weightparam(2) ...
                        + weightparam(3)*weights(tt-1) ...
                        + weightparam(4)*score_weights(tt-1);                
    else
        
        weights(tt) = (1-weightparam(3))*weightparam(2) ...
                        + weightparam(3)*weights(tt-1);                    
    end
% The commented code below is problematic, because it causes my_fmincon to 
% fail sometimes. If the true optimal likelihood lies at a point where the
% forecasts are non-pd, it tries to go closer and closer to this point,
% sometimes stopping at a point where there are non-finite derivatives,
% which causes my_fmincon to fail.
%     if any(eig(Sigma_(:,:,tt))<0) %This is if the forecasted Sigmas are not pd, which can happen despite all in-sample sigmas being pd. in this case likely the estimated scoreparam is too high for stationarity.
%         nLogL = inf;
%         logLcontr = NaN;
%         SigmaE = NaN;
%         varargout{1} = NaN;      
%         varargout{2} = NaN;
%         return
%     end
end
[Sigma_(:,:,T+1:T+t_ahead), pdadjustments_made] = makepd(Sigma_(:,:,T+1:T+t_ahead));
%% Log-Likelihood
nLogL = -sum(logLcontr);
weights_out.weights = weights;
weights_out.score_weights = score_weights;
%% Parameter Output
if nargout >= 5
    param_out.intrcpt = intrcpt;
    param_out.scoreparam = scoreparam;
    param_out.garch_param = garchparam;
    param_out.weight_param = weightparam;
    param_out.n_itRiesz2 = n_itRiesz2;
    param_out.nu_itRiesz2 = nu_itRiesz2;
    param_out.n_itRiesz2 = n_itRiesz2;
    param_out.nu_FRiesz = nu_FRiesz;
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
