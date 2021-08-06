function [ nLogL, logLcontr, SigmaE, ScaledScore1, varargout ] = ...
    gas_scalar_BEKK_alternative_scaling_likeRec( param, p, q, X, dist )
%
% Michael Stollenwerk
% michael.stollenwerk@live.com
% 11.08.2020

[k,~,T] = size(X);
k_ = k*(k+1)/2;
%% Parameters
intrcpt = ivechchol(param(1:k_));
scoreparam1 = param(k_ + 1 : k_ + p);
scoreparam2 = param(k_ + p + 1 : k_ + 2*p);
garchparam = param(k_ + 2*p : k_+ 2*p + q);

if strcmp( dist, 'Wish' )
    n = param(k_+ 2*p + q + 1);    
elseif strcmp( dist, 'iWish' )
    n = param(k_+ 2*p + q + 1);
elseif strcmp( dist, 'LaplaceWish' )
    n = param(k_+ 2*p + q + 1);      
elseif strcmp( dist, 'tWish' )
    n = param(k_+ 2*p + q + 1); 
    nu = param(k_+ 2*p + q + 2);
elseif strcmp( dist, 'itWish' )
    n = param(k_+ 2*p + q + 1); 
    nu = param(k_+ 2*p + q + 2);
elseif strcmp( dist, 'F' )
    n = param(k_+ 2*p + q + 1); 
    nu = param(k_+ 2*p + q + 2);
elseif strcmp( dist, 'Riesz' )
    n = param(k_+ 2*p + q + 1 : k_+ 2*p + q + k);    
elseif strcmp( dist, 'iRiesz' )
    n = param(k_+ 2*p + q + 1 : k_+ 2*p + q + k);  
elseif strcmp( dist, 'tRiesz' )
    n = param(k_+ 2*p + q + 1 : k_+ 2*p + q + k); 
    nu = param(k_+ p + q + k + 1);
elseif strcmp( dist, 'itRiesz' )
    n = param(k_+ 2*p + q + 1 : k_+ 2*p + q + k); 
    nu = param(k_+ 2*p + q + k + 1);
elseif strcmp( dist, 'FRiesz' )
    n = param(k_+ 2*p + q + 1 : k_+ 2*p + q + k); 
    nu = param(k_+ 2*p + q + k + 1 : k_+ 2*p + q + k + k);        
end

if nargout >= 5
    param_out.intrcpt = intrcpt;
    param_out.scoreparam1 = scoreparam1;
    param_out.scoreparam2 = scoreparam2;
    param_out.garch_param = garchparam;
    param_out.n = n;
    
    if exist('nu','var')
        param_out.nu = nu;
    end

    param_out.all = param;
    
    varargout{1} = param_out;
end
%% Data Storage
SigmaE = NaN(k,k,T+22);
ScaledScore1 = NaN(k,k,T);
ScaledScore2 = NaN(k,k,T);
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
G = Dmatrix(k);
iG = (G'*G)\G';
L = ELmatrix(k);
I = eye(k);
for tt=1:T
    SigmaE(:,:,tt) = intrcpt;
    for jj = 1:p
        if (tt-jj) > 0
            SigmaE(:,:,tt) = SigmaE(:,:,tt) ...
                                + scoreparam1(jj)*ScaledScore1(:,:,tt-jj) ...
                                + scoreparam2(jj)*ScaledScore2(:,:,tt-jj);
        end
    end
    for jj = 1:q
        if (tt-jj) <= 0
            SigmaE(:,:,tt) = SigmaE(:,:,tt) + garchparam(jj)*ini;
        else
            SigmaE(:,:,tt) = SigmaE(:,:,tt) + garchparam(jj)*SigmaE(:,:,tt-jj);
        end
    end
    
    try
        % Centering distributions
        if strcmp( dist, 'Wish' )
            Y = n;
            Sigma_ = SigmaE(:,:,tt)/Y;
            dSigmaEdSigma = Y;
        elseif strcmp( dist, 'LaplaceWish' )
            Y = n;
            Sigma_ = SigmaE(:,:,tt)/Y;
            dSigmaEdSigma = Y;
        elseif strcmp( dist, 'iWish' )
            Y = 1/(n-k-1);
            Sigma_ = SigmaE(:,:,tt)/Y;
            dSigmaEdSigma = Y;
        elseif strcmp( dist, 'tWish')
            Y = n*nu/(nu-2);
            Sigma_ = SigmaE(:,:,tt)/Y;
            dSigmaEdSigma = Y;
        elseif strcmp( dist, 'itWish')
            Y = 1/(n-k-1);
            Sigma_ = SigmaE(:,:,tt)/Y;
            dSigmaEdSigma = Y;
        elseif strcmp( dist, 'F' )
            Y = n/(nu-k-1);
            Sigma_ = SigmaE(:,:,tt)/Y;
            dSigmaEdSigma = Y;
        elseif strcmp( dist, 'Riesz' )
            Y = diag(n);
            CsigE = chol(SigmaE(:,:,tt),'lower');
            C = CsigE/sqrtm(Y);
            Sigma_ = C*C';
            dSigmaEdSigma = G'*kron(C*Y,I)*L'/(G'*kron(C,I)*L')*(G'*G);
        elseif strcmp( dist, 'iRiesz' )
            Y = inv(diag(n-k-1));
            iCdotSigE = inv(chol(inv(SigmaE(:,:,tt)),'lower'));
            iCdot = iCdotSigE/sqrtm(Y);
            Sigma_ = iCdot'*iCdot;
            dSigmaEdSigma = G'*kron(iCdot',iCdot'*Y*iCdot)*L'/(G'*kron(iCdot',Sigma_)*L')*(G'*G);
        elseif strcmp( dist, 'tRiesz' )
            Y = diag(n).*nu./(nu - 2);
            CsigE = chol(SigmaE(:,:,tt),'lower');
            C = CsigE/sqrtm(Y);
            Sigma_ = C*C';
            dSigmaEdSigma = G'*kron(C*Y,I)*L'/(G'*kron(C,I)*L')*(G'*G);
        elseif strcmp( dist, 'itRiesz' )
            Y = inv(diag(n-k-1));
            iCdotSigE = inv(chol(inv(SigmaE(:,:,tt)),'lower'));
            iCdot = iCdotSigE/sqrtm(Y);
            Sigma_ = iCdot'*iCdot;
            dSigmaEdSigma = G'*kron(iCdot',iCdot'*Y*iCdot)*L'/(G'*kron(iCdot',Sigma_)*L')*(G'*G);
        elseif strcmp( dist, 'FRiesz' )
            Y = matvFRieszexpmat(n,nu);
            CsigE = chol(SigmaE(:,:,tt),'lower');
            C = CsigE/sqrtm(Y);
            Sigma_ = C*C';
            dSigmaEdSigma = G'*kron(C*Y,I)*L'/(G'*kron(C,I)*L')*(G'*G);
        end
        % Likelihood Evaluation    
        if exist('nu','var')
            [logLcontr(tt), score, ~, ~, ~] = logpdf( dist, X(:,:,tt), Sigma_, n, nu );                           %%%%%%%%%%%%%%%%%%%%%%%
        else
            [logLcontr(tt), score, ~, ~, ~] = logpdf( dist, X(:,:,tt), Sigma_, n );
        end
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
    ScaledScore1(:,:,tt) = ivech(iG*kron2(Sigma_)*iG'/dSigmaEdSigma*score.Sigma_');
    ScaledScore2(:,:,tt) = ivech(iG*vec2(Sigma_)*iG'/dSigmaEdSigma*score.Sigma_');
    
end
%% Fcst
for tt=T+1:T+22
    SigmaE(:,:,tt) = intrcpt;
    for jj = 1:p
        if (tt-jj) <= T
            SigmaE(:,:,tt) = SigmaE(:,:,tt) ...
                                + scoreparam1(jj)*ScaledScore1(:,:,tt-jj) ...
                                + scoreparam2(jj)*ScaledScore2(:,:,tt-jj);
        end
    end
    for jj = 1:q
        SigmaE(:,:,tt) = SigmaE(:,:,tt) + garchparam(jj)*SigmaE(:,:,tt-jj);
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
