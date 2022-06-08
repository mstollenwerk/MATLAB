function [ Sigma_, S, pdadjustments_made ] = ...
    gas_scalar_BEKK_Rec( param, p, q, R, dist )
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
% ini = intrcpt./(1 - sum(garchparam));
%% Data Storage
Sigma_ = NaN(k,k,T+t_ahead);
S = NaN(k,k,T,2);
logLcontr = NaN(T,1);
%% Recursion
pdadjustments_made = NaN(T+t_ahead,1);
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
    
    [Sigma_(:,:,tt), pdadjustments_made(tt)] = makepd(Sigma_(:,:,tt));
    % Likelihood Evaluation         
    [~, logLcontr(tt), score] = ...
        matvsLogLike( dist, Sigma_(:,:,tt), param(k_+2*p+q+1:end), R(:,:,tt) );   
        
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
    [Sigma_(:,:,tt), pdadjustments_made(tt)] = makepd(Sigma_(:,:,tt));
end

