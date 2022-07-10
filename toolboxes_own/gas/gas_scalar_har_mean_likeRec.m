function [ nLogL, logLcontr, dyn, S, varargout ] = ...
    gas_scalar_har_mean_likeRec( param, R, dist )
%
% Michael Stollenwerk
% michael.stollenwerk@live.com
% 19.06.2022

t_ahead = 220;

[p,~,T] = size(R);
p_ = p*(p+1)/2;

par = gas_scalar_har_mean_param_transform(dist,param,p);

s1Sig = par.Sig.score1;
s2Sig = par.Sig.score2;
gDSig = par.Sig.garchD;
gWSig = par.Sig.garchW;
gMSig = par.Sig.garchM;
theta = param(p_+6:end);

% Initialize recursion at unconditional mean (stationarity assumed).
% iniMatSig = diag(1./sqrt(1-par.Sig.garchD.^2-par.Sig.garchW.^2-par.Sig.garchM.^2));
% iniSig = iniMatSig*par.Sig.intrcpt*iniMatSig;
% Initialize with Backcast. Inspired by Sheppard MFE Toolbox.
m = ceil(sqrt(T));
w = .06 * .94.^(0:(m-1));
w = reshape(w/sum(w),[1 1 m]);
backCast = sum(bsxfun(@times,w, R(:,:,1:m)),3);
iniSig = backCast;

Sig = NaN(p,p,T+t_ahead);
Sig_score = NaN(p,p,T,2);

logLcontr = NaN(T,1);
%% Recursion
for tt=1:T
    
    % Recursion                 
    Sig(:,:,tt) = par.Sig.intrcpt;    
    for jj = 1:22
        if (tt-jj)<1
            if jj == 1
                Sig(:,:,tt) = Sig(:,:,tt) ...
                                 + gDSig*iniSig ...
                                 + 1/5*gWSig*iniSig ...
                                 + 1/22*gMSig*iniSig;
            elseif (1 < jj) && (jj <= 5)
                Sig(:,:,tt) = Sig(:,:,tt) ...
                                 + 1/5*gWSig*iniSig ...
                                 + 1/22*gMSig*iniSig; 
            elseif (5 < jj) && (jj <= 22)
                Sig(:,:,tt) = Sig(:,:,tt) ...
                                 + 1/22*gMSig*iniSig;
            end  
        else
            if jj == 1
                Sig(:,:,tt) = Sig(:,:,tt) ...
                                 + s1Sig*Sig_score(:,:,tt-jj,1) ...
                                 + s2Sig*Sig_score(:,:,tt-jj,2) ...
                                 + gDSig*Sig(:,:,tt-jj) ...
                                 + 1/5*gWSig*Sig(:,:,tt-jj) ...
                                 + 1/22*gMSig*Sig(:,:,tt-jj);
            elseif (1 < jj) && (jj <= 5)
                Sig(:,:,tt) = Sig(:,:,tt) ...
                                 + 1/5*gWSig*Sig(:,:,tt-jj) ...
                                 + 1/22*gMSig*Sig(:,:,tt-jj);
            elseif (5 < jj) && (jj <= 22)
                Sig(:,:,tt) = Sig(:,:,tt) ...
                                 + 1/22*gMSig*Sig(:,:,tt-jj);
            end
        end
    end
    
    try
        % Likelihood Evaluation         
        [~, logLcontr(tt), score ] = ...
            matvsLogLike( dist, Sig(:,:,tt), theta, R(:,:,tt) );   
    catch ME
%         tt
%         ME.message
%         [ME.stack.line]'
%         {ME.stack.name}'
        nLogL = inf;
        logLcontr = NaN;
        dyn = NaN;
        S = NaN;
        varargout{1} = NaN;      
        varargout{2} = NaN;
        return
    end
    
    % Scaled Scores 
    SigDel = Sig(:,:,tt)*score.SigmaNonSym;
    SigDelSig = SigDel*Sig(:,:,tt);
    Sig_score(:,:,tt,1) = SigDelSig+SigDelSig';
    Sig_score(:,:,tt,2) = trace(SigDel)*Sig(:,:,tt);
    
end
%% Fcst

for tt=T+1
    Sig(:,:,tt) = par.Sig.intrcpt;
    for jj = 1:22
        if jj == 1
            Sig(:,:,tt) = Sig(:,:,tt) ...
                             + s1Sig*Sig_score(:,:,tt-jj,1) ...
                             + s2Sig*Sig_score(:,:,tt-jj,2) ...
                             + gDSig*Sig(:,:,tt-jj) ...
                             + 1/5*gWSig*Sig(:,:,tt-jj) ...
                             + 1/22*gMSig*Sig(:,:,tt-jj);
        elseif (1 < jj) && (jj <= 5)
            Sig(:,:,tt) = Sig(:,:,tt) ...
                             + 1/5*gWSig*Sig(:,:,tt-jj) ...
                             + 1/22*gMSig*Sig(:,:,tt-jj);
        elseif (5 < jj) && (jj <= 22)
            Sig(:,:,tt) = Sig(:,:,tt) ...
                             + 1/22*gMSig*Sig(:,:,tt-jj);
        end
    end
end

for tt=T+2:T+t_ahead
    Sig(:,:,tt) = par.Sig.intrcpt;
    for jj = 1:22    
        if jj == 1
            Sig(:,:,tt) = Sig(:,:,tt) ...
                             + gDSig*Sig(:,:,tt-jj) ...
                             + 1/5*gWSig*Sig(:,:,tt-jj) ...
                             + 1/22*gMSig*Sig(:,:,tt-jj);
        elseif (1 < jj) && (jj <= 5)
            Sig(:,:,tt) = Sig(:,:,tt) ...
                             + 1/5*gWSig*Sig(:,:,tt-jj) ...
                             + 1/22*gMSig*Sig(:,:,tt-jj);
        elseif (5 < jj) && (jj <= 22)
            Sig(:,:,tt) = Sig(:,:,tt) ...
                             + 1/22*gMSig*Sig(:,:,tt-jj);
        end
    end
end

    
% The commented code below is problematic, because it causes my_fmincon to 
% fail sometimes. If the true optimal likelihood lies at a point where the
% forecasts are non-pd, it tries to go closer and closer to this point,
% sometimes stopping at a point where there are non-finite derivatives,
% which causes my_fmincon to fail.
%     if any(eig(Sig(:,:,tt))<0) %This is if the forecasted Sigmas are not pd, which can happen despite all in-sample sigmas being pd. in this case likely the estimated scoreparam is too high for stationarity.
%         nLogL = inf;
%         logLcontr = NaN;
%         Sig = NaN;
%         varargout{1} = NaN;      
%         varargout{2} = NaN;
%         return
%     end
% end
% [Sigma_(:,:,T+1:T+t_ahead), pdadjustments_made] = makepd(Sigma_(:,:,T+1:T+t_ahead));
%% Output
nLogL = -sum(logLcontr);
dyn.Sig = Sig;
S.Sig = Sig_score;
%% Parameter Output
if nargout >= 5
    param_out = par;
    param_out.all = param;
    
    varargout{1} = param_out;
end

end
