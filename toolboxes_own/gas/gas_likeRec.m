function [ nLogL, logLcontr, dyn, S, varargout ] = ...
    gas_likeRec( param, R, dist )
%
% Michael Stollenwerk
% michael.stollenwerk@live.com
% 19.06.2022

t_ahead = 220;

[p,~,T] = size(R);
p_ = p*(p+1)/2;

par = gas_param_transform(dist,param,p);
sSig1 = diag(par.Sig.score1);
sSig2 = diag(par.Sig.score2);
gDSig = diag(par.Sig.garchD);
gWSig = diag(par.Sig.garchW);
gMSig = diag(par.Sig.garchM);

existN = isfield(par,'n');
existNu = isfield(par,'nu');

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
Ssig = NaN(p,p,T,2);

if existN
    n = NaN(length(par.n.intrcpt),T+t_ahead);
    Sn = NaN(length(par.n.intrcpt),T);
    iniN = par.n.intrcpt ./ (1-par.n.garchD-par.n.garchW-par.n.garchM);
end
if existNu
    nu = NaN(length(par.nu.intrcpt),T+t_ahead);
    Snu = NaN(length(par.nu.intrcpt),T);
    iniNu = par.nu.intrcpt ./ (1-par.nu.garchD-par.nu.garchW-par.nu.garchM);
end

logLcontr = NaN(T,1);
%% Recursion
for tt=1:T
    
    % Recursion                 
    Sig(:,:,tt) = par.Sig.intrcpt;
    if existN
        n(:,tt) = par.n.intrcpt;
    end
    if existNu
        nu(:,tt) = par.nu.intrcpt;
    end
    for jj = 1:22
        if (tt-jj)<1
            if jj == 1
                Sig(:,:,tt) = Sig(:,:,tt) ...
                                 + gDSig*iniSig*gDSig ...
                                 + 1/5*gWSig*iniSig*gWSig ...
                                 + 1/22*gMSig*iniSig*gMSig;
                if existN
                    n(:,tt) = n(:,tt) ...
                              + par.n.garchD*iniN ...
                              + 1/5*par.n.garchW*iniN ...
                              + 1/22*par.n.garchM*iniN;
                end
                if existNu
                    nu(:,tt) = nu(:,tt) ...
                              + par.nu.garchD*iniNu ...
                              + 1/5*par.nu.garchW*iniNu ...
                              + 1/22*par.nu.garchM*iniNu;
                end
            elseif (1 < jj) && (jj <= 5)
                Sig(:,:,tt) = Sig(:,:,tt) ...
                                 + 1/5*gWSig*iniSig*gWSig ...
                                 + 1/22*gMSig*iniSig*gMSig; 
                if existN
                    n(:,tt) = n(:,tt) ...
                              + 1/5*par.n.garchW*iniN ...
                              + 1/22*par.n.garchM*iniN;
                end
                if existNu
                    nu(:,tt) = nu(:,tt) ...
                              + 1/5*par.nu.garchW*iniNu ...
                              + 1/22*par.nu.garchM*iniNu;
                end
            elseif (5 < jj) && (jj <= 22)
                Sig(:,:,tt) = Sig(:,:,tt) ...
                                 + 1/22*gMSig*iniSig*gMSig;  
                if existN
                    n(:,tt) = n(:,tt) ...
                              + 1/22*par.n.garchM*iniN;
                end
                if existNu
                    nu(:,tt) = nu(:,tt) ...
                              + 1/22*par.nu.garchM*iniNu;
                end
            end           
        else
            if jj == 1
                Sig(:,:,tt) = Sig(:,:,tt) ...
                                 + sSig1*Ssig(:,:,tt-jj,1)*sSig1 ...
                                 + sSig2*Ssig(:,:,tt-jj,2)*sSig2 ...
                                 + gDSig*Sig(:,:,tt-jj)*gDSig ...
                                 + 1/5*gWSig*Sig(:,:,tt-jj)*gWSig ...
                                 + 1/22*gMSig*Sig(:,:,tt-jj)*gMSig;
                if existN
                    n(:,tt) = n(:,tt) ...
                              + par.n.score*Sn(:,tt-jj) ...
                              + par.n.garchD*n(:,tt-jj) ...
                              + 1/5*par.n.garchW*n(:,tt-jj) ...
                              + 1/22*par.n.garchM*n(:,tt-jj);
                end
                if existNu
                    nu(:,tt) = nu(:,tt) ...
                              + par.nu.score*Snu(:,tt-jj) ...
                              + par.nu.garchD*nu(:,tt-jj) ...
                              + 1/5*par.nu.garchW*nu(:,tt-jj) ...
                              + 1/22*par.nu.garchM*nu(:,tt-jj);
                end
            elseif (1 < jj) && (jj <= 5)
                Sig(:,:,tt) = Sig(:,:,tt) ...
                                 + 1/5*gWSig*Sig(:,:,tt-jj)*gWSig ...
                                 + 1/22*gMSig*Sig(:,:,tt-jj)*gMSig;                
                if existN
                    n(:,tt) = n(:,tt) ...
                              + 1/5*par.n.garchW*n(:,tt-jj) ...
                              + 1/22*par.n.garchM*n(:,tt-jj);
                end
                if existNu
                    nu(:,tt) = nu(:,tt) ...
                              + 1/5*par.nu.garchW*nu(:,tt-jj) ...
                              + 1/22*par.nu.garchM*nu(:,tt-jj);
                end
            elseif (5 < jj) && (jj <= 22)
                Sig(:,:,tt) = Sig(:,:,tt) ...
                                 + 1/22*gMSig*Sig(:,:,tt-jj)*gMSig;                   
                if existN
                    n(:,tt) = n(:,tt) ...
                              + 1/22*par.n.garchM*n(:,tt-jj);
                end
                if existNu
                    nu(:,tt) = nu(:,tt) ...
                              + 1/22*par.nu.garchM*nu(:,tt-jj);
                end
            end
        end
    end
    
    if existN && existNu
        theta = [n(:,tt); nu(:,tt)];
    elseif existN 
        theta = n(:,tt);
    elseif existNu
        theta = nu(:,tt);
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
    Ssig(:,:,tt,1) = SigDelSig+SigDelSig';
    Ssig(:,:,tt,2) = trace(SigDel)*Sig(:,:,tt);
    
    if existN
        Sn(:,tt) = score.n_originalpdf_scaled;
    end
    if existNu
        Snu(:,tt) = score.nu_originalpdf_scaled;
    end
    
end
%% Fcst

for tt=T+1
    Sig(:,:,tt) = par.Sig.intrcpt;
    if existN
        n(:,tt) = par.n.intrcpt;
    end
    if existNu
        nu(:,tt) = par.nu.intrcpt;
    end
    for jj = 1:22
        if jj == 1
            Sig(:,:,tt) = Sig(:,:,tt) ...
                             + sSig1*Ssig(:,:,tt-jj,1)*sSig1 ...
                             + sSig2*Ssig(:,:,tt-jj,2)*sSig2 ...
                             + gDSig*Sig(:,:,tt-jj)*gDSig ...
                             + 1/5*gWSig*Sig(:,:,tt-jj)*gWSig ...
                             + 1/22*gMSig*Sig(:,:,tt-jj)*gMSig;
            if existN
                n(:,tt) = n(:,tt) ...
                          + par.n.score*Sn(:,tt-jj) ...
                          + par.n.garchD*n(:,tt-jj) ...
                          + 1/5*par.n.garchW*n(:,tt-jj) ...
                          + 1/22*par.n.garchM*n(:,tt-jj);
            end
            if existNu
                nu(:,tt) = nu(:,tt) ...
                          + par.nu.score*Snu(:,tt-jj) ...
                          + par.nu.garchD*nu(:,tt-jj) ...
                          + 1/5*par.nu.garchW*nu(:,tt-jj) ...
                          + 1/22*par.nu.garchM*nu(:,tt-jj);
            end
        elseif (1 < jj) && (jj <= 5)
            Sig(:,:,tt) = Sig(:,:,tt) ...
                             + 1/5*gWSig*Sig(:,:,tt-jj)*gWSig ...
                             + 1/22*gMSig*Sig(:,:,tt-jj)*gMSig;                
            if existN
                n(:,tt) = n(:,tt) ...
                          + 1/5*par.n.garchW*n(:,tt-jj) ...
                          + 1/22*par.n.garchM*n(:,tt-jj);
            end
            if existNu
                nu(:,tt) = nu(:,tt) ...
                          + 1/5*par.nu.garchW*nu(:,tt-jj) ...
                          + 1/22*par.nu.garchM*nu(:,tt-jj);
            end
        elseif (5 < jj) && (jj <= 22)
            Sig(:,:,tt) = Sig(:,:,tt) ...
                             + 1/22*gMSig*Sig(:,:,tt-jj)*gMSig;                   
            if existN
                n(:,tt) = n(:,tt) ...
                          + 1/22*par.n.garchM*n(:,tt-jj);
            end
            if existNu
                nu(:,tt) = nu(:,tt) ...
                          + 1/22*par.nu.garchM*nu(:,tt-jj);
            end
        end
    end
end

for tt=T+2:T+t_ahead
    Sig(:,:,tt) = par.Sig.intrcpt;
    if existN
        n(:,tt) = par.n.intrcpt;
    end
    if existNu
        nu(:,tt) = par.nu.intrcpt;
    end
    for jj = 1:22    
        if jj == 1
            Sig(:,:,tt) = Sig(:,:,tt) ...
                             + gDSig*Sig(:,:,tt-jj)*gDSig ...
                             + 1/5*gWSig*Sig(:,:,tt-jj)*gWSig ...
                             + 1/22*gMSig*Sig(:,:,tt-jj)*gMSig;
            if existN
                n(:,tt) = n(:,tt) ...
                          + par.n.garchD*n(:,tt-jj) ...
                          + 1/5*par.n.garchW*n(:,tt-jj) ...
                          + 1/22*par.n.garchM*n(:,tt-jj);
            end
            if existNu
                nu(:,tt) = nu(:,tt) ...
                          + par.nu.garchD*nu(:,tt-jj) ...
                          + 1/5*par.nu.garchW*nu(:,tt-jj) ...
                          + 1/22*par.nu.garchM*nu(:,tt-jj);
            end
        elseif (1 < jj) && (jj <= 5)
            Sig(:,:,tt) = Sig(:,:,tt) ...
                             + 1/5*gWSig*Sig(:,:,tt-jj)*gWSig ...
                             + 1/22*gMSig*Sig(:,:,tt-jj)*gMSig;                
            if existN
                n(:,tt) = n(:,tt) ...
                          + 1/5*par.n.garchW*n(:,tt-jj) ...
                          + 1/22*par.n.garchM*n(:,tt-jj);
            end
            if existNu
                nu(:,tt) = nu(:,tt) ...
                          + 1/5*par.nu.garchW*nu(:,tt-jj) ...
                          + 1/22*par.nu.garchM*nu(:,tt-jj);
            end
        elseif (5 < jj) && (jj <= 22)
            Sig(:,:,tt) = Sig(:,:,tt) ...
                             + 1/22*gMSig*Sig(:,:,tt-jj)*gMSig;                   
            if existN
                n(:,tt) = n(:,tt) ...
                          + 1/22*par.n.garchM*n(:,tt-jj);
            end
            if existNu
                nu(:,tt) = nu(:,tt) ...
                          + 1/22*par.nu.garchM*nu(:,tt-jj);
            end
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
S.Sig = Ssig;
if existNu
    dyn.nu = nu;
    S.nu = Snu;
end
if existN
    dyn.n = n;
    S.n = Sn;
end
%% Parameter Output
if nargout >= 5
    param_out = par;
    param_out.all = param;
    
    varargout{1} = param_out;
end

end
