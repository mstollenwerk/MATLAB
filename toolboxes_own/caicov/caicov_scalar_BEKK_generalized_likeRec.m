function [ nLogL, logLcontr, SigmaE, weights_out, varargout ] = ...
    caicov_scalar_BEKK_generalized_likeRec( param, p, q, X )
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
% weightparam = param(k_ + p + q + 1 : k_ + p + q + 5);
% n_itRiesz2 = param(k_ + p + q + 6 : k_ + p + q + 5 + k);
% nu_itRiesz2 = param(k_ + p + q + 5 + k + 1);
% n_FRiesz = param(k_ + p + q + 5 + k + 2 : k_ + p + q + 5 + k + 1 + k);
% nu_FRiesz = param(k_ + p + q + 5 + k + 1 + k + 1: k_ + p + q + 5 + k + 1 + k + k);
weightparam = param(k_ + p + q + 1 : k_ + p + q + 4);
n_itRiesz2 = param(k_ + p + q + 5 : k_ + p + q + 4 + k);
nu_itRiesz2 = param(k_ + p + q + 4 + k + 1);
n_FRiesz = param(k_ + p + q + 4 + k + 2 : k_ + p + q + 4 + k + 1 + k);
nu_FRiesz = param(k_ + p + q + 4 + k + 1 + k + 1: k_ + p + q + 4 + k + 1 + k + k);

if nargout >= 4
    param_out.intrcpt = intrcpt;
    param_out.archparam = archparam;
    param_out.garchparam = garchparam;

    param_out.all = param;
    
    varargout{1} = param_out;
end
%% Data Storage
SigmaE = NaN(k,k,T+t_ahead);
weights = NaN(T+t_ahead,1);
logLcontr = NaN(T,1);
Lcontr_itRiesz2 = NaN(T,1);
Lcontr_FRiesz = NaN(T,1);
score_weights = NaN(T+t_ahead,1);
%% Recursion
weights(1) = weightparam(1);
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
    if tt>1
%         [~, relnegcovs_tt1, ~] = negcovs(X(:,:,tt-1));
%         ldS = logdet(SigmaE(:,:,tt));
%         if ldS < 0 
%             weights(tt) = weightparam(2) ...  
%                             + weightparam(3)*weights(tt-1) ...
%                             + weightparam(4)*ldS;
%     %                         + weightparam(5)*sum(diag(SigmaE(:,:,tt)))^(1/20);
%         else
%             weights(tt) = weightparam(2) ...  
%                             + weightparam(3)*weights(tt-1) ...
%                             + weightparam(5)*ldS;
%         end
%         weights(tt) = weightparam(2) ...
%                         + weightparam(3)*weights(tt-1) ...
%                         + weightparam(4)*log(det(SigmaE(:,:,tt))+1);
%                     
        weights(tt) = (1-weightparam(3))*weightparam(2) ...
                        + weightparam(3)*weights(tt-1) ...
                        + weightparam(4)*score_weights(tt-1);
    end   

    try
        % Centering distributions
        Sigma_itRiesz2 = matvStandardize('itRiesz2', SigmaE(:,:,tt), [n_itRiesz2; nu_itRiesz2]);
        Sigma_FRiesz = matvStandardize('FRiesz', SigmaE(:,:,tt), [n_FRiesz; nu_FRiesz]);
        % Likelihood Evaluation
        nLogLcontr_itRiesz2 = matvitRiesz2like(Sigma_itRiesz2,n_itRiesz2,nu_itRiesz2,X(:,:,tt));
        nLogLcontr_FRiesz = matvFRieszlike(Sigma_FRiesz,n_FRiesz,nu_FRiesz,X(:,:,tt));
        Lcontr_itRiesz2(tt) = exp(-nLogLcontr_itRiesz2);
        Lcontr_FRiesz(tt) = exp(-nLogLcontr_FRiesz);
        Like_tt = weights(tt)*Lcontr_itRiesz2(tt)+ (1-weights(tt))*Lcontr_FRiesz(tt);
        score_weights(tt) = (Lcontr_itRiesz2(tt) - Lcontr_FRiesz(tt))/Like_tt;
        logLcontr(tt) = log(Like_tt);
        
%         if weights(tt) > 1 || weights(tt) < 0
%             nLogL = inf;
%             logLcontr = NaN;
%             SigmaE = NaN;
%             weights_out = NaN;
%             varargout{1} = NaN;      
%             varargout{2} = NaN;
%             return
%         end
            
        
    catch ME
%         tt
%         ME.message
%         [ME.stack.line]'
%         {ME.stack.name}'
        nLogL = inf;
        logLcontr = NaN;
        SigmaE = NaN;
        weights_out = NaN;        
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
    if tt == T+1
%         ldS = logdet(SigmaE(:,:,tt));
%         if ldS < 0 
%             weights(tt) = weightparam(2) ...  
%                             + weightparam(3)*weights(tt-1) ...
%                             + weightparam(4)*ldS;
% %                             + weightparam(5)*sum(diag(SigmaE(:,:,tt)))^(1/20);
%         else
%             weights(tt) = weightparam(2) ...  
%                             + weightparam(3)*weights(tt-1) ...
%                             + weightparam(5)*ldS;
%         end

        weights(tt) = (1-weightparam(3))*weightparam(2) ...
                        + weightparam(3)*weights(tt-1) ...
                        + weightparam(4)*score_weights(tt-1);
                    
%         weights(tt) = weightparam(2) ...
%                         + weightparam(3)*weights(tt-1) ...
%                         + weightparam(4)*log(det(SigmaE(:,:,tt))+1);                    
    else
        
        weights(tt) = (1-weightparam(3))*weightparam(2) ...
                        + weightparam(3)*weights(tt-1); 
                    
%         weights(tt) = weightparam(2) ... %Same ase for tt== T+1, did in development, can be changed.
%                         + weightparam(3)*weights(tt-1) ...
%                         + weightparam(4)*log(det(SigmaE(:,:,tt))+1);                     
    end
end
%% Log-Likelihood
nLogL = -sum(logLcontr);
weights_out.weights = weights;
weights_out.score_weights = score_weights;
%% Fit-Plot(s)
if nargout >= 6
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
