function [ nLogL, logLcontr, dyn, S, varargout ] = ...
    gas_scalar_likeRec_( param, R, dist )
%
% Michael Stollenwerk
% michael.stollenwerk@live.com
% 29.11.2022

t_ahead = 220;
    
[p,~,T] = size(R);
p_ = p*(p+1)/2;

par = gas_scalar_param_transform_(dist,param,p);

s1Sig = par.Sig.score1;
s2Sig = par.Sig.score2;
gDSig = par.Sig.garchD;
gWSig = par.Sig.garchW;
gMSig = par.Sig.garchM;

existN = isfield(par,'n');
existNu = isfield(par,'nu');

if existN
    n = NaN(length(par.n.intrcpt),T+t_ahead);
    n_score = NaN(length(par.n.intrcpt),T);
    s_n = par.n.score;
    gDn = par.n.garchD;
    gWn = par.n.garchW;
    gMn = par.n.garchM;
    i_n = par.n.intrcpt;
    ini_n = i_n/(1-gDn-gWn-gMn);
end
if existNu
    nu = NaN(length(par.nu.intrcpt),T+t_ahead);
    nu_score = NaN(length(par.nu.intrcpt),T);
    s_nu = par.nu.score;
    gDnu = par.nu.garchD;
    gWnu = par.nu.garchW;
    gMnu = par.nu.garchM;
    i_nu = par.nu.intrcpt;
    ini_nu = i_nu/(1-gDnu-gWnu-gMnu);
end

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
    
    % Recursion Sigma                 
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
    
    % Recursion n 
    if existN          
        n(:,tt) = i_n;    
        for jj = 1:22
            if (tt-jj)<1
                if jj == 1
                    n(:,tt) = n(:,tt) ...
                                     + gDn*ini_n ...
                                     + 1/5*gWn*ini_n ...
                                     + 1/22*gMn*ini_n;
                elseif (1 < jj) && (jj <= 5)
                    n(:,tt) = n(:,tt) ...
                                     + 1/5*gWn*ini_n ...
                                     + 1/22*gMn*ini_n; 
                elseif (5 < jj) && (jj <= 22)
                    n(:,tt) = n(:,tt) ...
                                     + 1/22*gMn*ini_n;
                end  
            else
                if jj == 1
                    n(:,tt) = n(:,tt) ...
                                     + s_n*n_score(:,tt-jj) ...
                                     + gDn*n(:,tt-jj) ...
                                     + 1/5*gWn*n(:,tt-jj) ...
                                     + 1/22*gMn*n(:,tt-jj);
                elseif (1 < jj) && (jj <= 5)
                    n(:,tt) = n(:,tt) ...
                                     + 1/5*gWn*n(:,tt-jj) ...
                                     + 1/22*gMn*n(:,tt-jj);
                elseif (5 < jj) && (jj <= 22)
                    n(:,tt) = n(:,tt) ...
                                     + 1/22*gMn*n(:,tt-jj);
                end
            end
        end
    end
    % Recursion nu 
    if existNu          
        nu(:,tt) = i_nu;    
        for jj = 1:22
            if (tt-jj)<1
                if jj == 1
                    nu(:,tt) = nu(:,tt) ...
                                     + gDnu*ini_nu ...
                                     + 1/5*gWnu*ini_nu ...
                                     + 1/22*gMnu*ini_nu;
                elseif (1 < jj) && (jj <= 5)
                    nu(:,tt) = nu(:,tt) ...
                                     + 1/5*gWnu*ini_nu ...
                                     + 1/22*gMnu*ini_nu; 
                elseif (5 < jj) && (jj <= 22)
                    nu(:,tt) = nu(:,tt) ...
                                     + 1/22*gMnu*ini_nu;
                end  
            else
                if jj == 1
                    nu(:,tt) = nu(:,tt) ...
                                     + s_nu*nu_score(:,tt-jj) ...
                                     + gDnu*nu(:,tt-jj) ...
                                     + 1/5*gWnu*nu(:,tt-jj) ...
                                     + 1/22*gMnu*nu(:,tt-jj);
                elseif (1 < jj) && (jj <= 5)
                    nu(:,tt) = nu(:,tt) ...
                                     + 1/5*gWnu*nu(:,tt-jj) ...
                                     + 1/22*gMnu*nu(:,tt-jj);
                elseif (5 < jj) && (jj <= 22)
                    nu(:,tt) = nu(:,tt) ...
                                     + 1/22*gMnu*nu(:,tt-jj);
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
    Sig_score(:,:,tt,1) = SigDelSig+SigDelSig';
    Sig_score(:,:,tt,2) = trace(SigDel)*Sig(:,:,tt);
    
    if existN
        n_score(:,tt) = score.n_originalpdf_scaled;
    end
    if existNu
        nu_score(:,tt) = score.nu_originalpdf_scaled;
    end
    
end
%% Fcst

for tt=T+1
    % Fcst Sigma
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
    % Fcst n
    if existN          
        n(:,tt) = i_n;    
        for jj = 1:22
            if jj == 1
                n(:,tt) = n(:,tt) ...
                                 + s_n*n_score(:,tt-jj) ...
                                 + gDn*n(:,tt-jj) ...
                                 + 1/5*gWn*n(:,tt-jj) ...
                                 + 1/22*gMn*n(:,tt-jj);
            elseif (1 < jj) && (jj <= 5)
                n(:,tt) = n(:,tt) ...
                                 + 1/5*gWn*n(:,tt-jj) ...
                                 + 1/22*gMn*n(:,tt-jj);
            elseif (5 < jj) && (jj <= 22)
                n(:,tt) = n(:,tt) ...
                                 + 1/22*gMn*n(:,tt-jj);
            end
        end
    end
    % Fcst nu
    if existNu          
        nu(:,tt) = i_nu;    
        for jj = 1:22
            if jj == 1
                nu(:,tt) = nu(:,tt) ...
                                 + s_nu*nu_score(:,tt-jj) ...
                                 + gDnu*nu(:,tt-jj) ...
                                 + 1/5*gWnu*nu(:,tt-jj) ...
                                 + 1/22*gMnu*nu(:,tt-jj);
            elseif (1 < jj) && (jj <= 5)
                nu(:,tt) = nu(:,tt) ...
                                 + 1/5*gWnu*nu(:,tt-jj) ...
                                 + 1/22*gMnu*nu(:,tt-jj);
            elseif (5 < jj) && (jj <= 22)
                nu(:,tt) = nu(:,tt) ...
                                 + 1/22*gMnu*nu(:,tt-jj);
            end
        end
    end       
end

for tt=T+2:T+t_ahead
    % Fcst Sigma
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
    % Fcst n
    if existN          
        n(:,tt) = i_n;    
        for jj = 1:22
            if jj == 1
                n(:,tt) = n(:,tt) ...
                                 + gDn*n(:,tt-jj) ...
                                 + 1/5*gWn*n(:,tt-jj) ...
                                 + 1/22*gMn*n(:,tt-jj);
            elseif (1 < jj) && (jj <= 5)
                n(:,tt) = n(:,tt) ...
                                 + 1/5*gWn*n(:,tt-jj) ...
                                 + 1/22*gMn*n(:,tt-jj);
            elseif (5 < jj) && (jj <= 22)
                n(:,tt) = n(:,tt) ...
                                 + 1/22*gMn*n(:,tt-jj);
            end
        end
    end
    % Fcst nu
    if existNu          
        nu(:,tt) = i_nu;    
        for jj = 1:22
            if jj == 1
                nu(:,tt) = nu(:,tt) ...
                                 + gDnu*nu(:,tt-jj) ...
                                 + 1/5*gWnu*nu(:,tt-jj) ...
                                 + 1/22*gMnu*nu(:,tt-jj);
            elseif (1 < jj) && (jj <= 5)
                nu(:,tt) = nu(:,tt) ...
                                 + 1/5*gWnu*nu(:,tt-jj) ...
                                 + 1/22*gMnu*nu(:,tt-jj);
            elseif (5 < jj) && (jj <= 22)
                nu(:,tt) = nu(:,tt) ...
                                 + 1/22*gMnu*nu(:,tt-jj);
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
S.Sig = Sig_score;
if existNu
    dyn.nu = nu;
    S.nu = nu_score;
end
if existN
    dyn.n = n;
    S.n = n_score;
end
%% Parameter Output
if nargout >= 5
    param_out = par;
    param_out.all = param;
    
    varargout{1} = param_out;
end

end
