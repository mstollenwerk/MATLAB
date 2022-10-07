function [ nLogL, logLcontr, Sigma_, varargout ] = ...
    caicov_scalar_BEKK_likeRec( param, p, q, R, dist )
%
% Michael Stollenwerk
% michael.stollenwerk@live.com
% 05.04.2021

t_ahead = 220;

[k,~,T] = size(R);
k_ = k*(k+1)/2;
%% Parameters
% Initialize with Backcast. Code copied from Sheppard MFE Toolbox.
m = ceil(sqrt(T));
w = .06 * .94.^(0:(m-1));
w = reshape(w/sum(w),[1 1 m]);
backCast = sum(bsxfun(@times,w,R(:,:,1:m)),3);
% % Initialize recursion at unconditional mean (stationarity assumed).
% ini = intrcpt./(1 - sum(garchparam) - sum(archparam));

intrcpt = ivechchol(param(1:k_));
archparam = param(k_ + 1 : k_ + p);
garchparam = param(k_ + p + 1 : k_+ p + q);
theta =  param(k_+p+q+1:end);
M = matvExpMat(dist,theta,k);
sqrtM = sqrt(M);
if any(diag(M)<0) %Checking existence of mean.
    nLogL = inf;
    logLcontr = NaN;
    Sigma_ = NaN;
    varargout{1} = NaN; 
    return
end

if nargout >= 4
    param_out.intrcpt = intrcpt;
    param_out.archparam = archparam;
    param_out.garchparam = garchparam;

    param_out.all = param;
    
    varargout{1} = param_out;
end
%% Data Storage
Sigma_ = NaN(k,k,T+t_ahead);
Cz = NaN(k,k,T);
%% Recursion

for tt=1:T
    
    Sigma_(:,:,tt) = intrcpt;
    for jj = 1:p
        if (tt-jj) <= 0
            Sigma_(:,:,tt) = Sigma_(:,:,tt) + archparam(jj)*backCast;
        else
            Sigma_(:,:,tt) = Sigma_(:,:,tt) + archparam(jj)*R(:,:,tt-jj);
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
        C = chol(Sigma_(:,:,tt),'lower');
        Cr = chol(R(:,:,tt),'lower');
		COm = C/sqrtM;
        Cz(:,:,tt) = COm\Cr;
        if any(strcmp( dist, {'iWish','iRiesz2','itWish','itRiesz2'} ))
            X(tt,:) = sum(sum(inv(Cz(:,:,tt)).^2));
        elseif strcmp( dist, 'F' )
            X(tt,:) = logdet(eye(k) + Cz(:,:,tt)*Cz(:,:,tt)');
        elseif strcmp( dist, 'FRiesz' )
            X(tt,:) = diag(chol(eye(k) + Cz(:,:,tt)*Cz(:,:,tt)','lower'));
        elseif strcmp( dist, 'iFRiesz2' )
            C_Z_twisted_t = Cr'/COm';
            B = eye(k) + C_Z_twisted_t*C_Z_twisted_t';
            X(tt,:) = diag(cholU(B));            
% 			C_Z_twisted_t = Cr/COm;
% 			X(tt,:) = diag(chol(inv(eye(k) + C_Z_twisted_t'*C_Z_twisted_t),'lower'));
        else
            X = [];
        end        
    catch ME
    %         tt
    %         ME.message
    %         [ME.stack.line]'
    %         {ME.stack.name}'
        nLogL = inf;
        logLcontr = NaN;
        Sigma_ = NaN;
        varargout{1} = NaN; 
        return
    end
    
end
%% Likelihood Evaluation
[nLogL, logLcontr] = matvLogLikeCz( dist, Cz, theta', X );
%% Fcst
for tt=T+1:T+t_ahead
    Sigma_(:,:,tt) = intrcpt;
    for jj = 1:p
        if (tt-jj) > T
            Sigma_(:,:,tt) = Sigma_(:,:,tt) + archparam(jj)*Sigma_(:,:,tt-jj);
        else
            Sigma_(:,:,tt) = Sigma_(:,:,tt) + archparam(jj)*R(:,:,tt-jj);
        end
    end
    for jj = 1:q
        Sigma_(:,:,tt) = Sigma_(:,:,tt) + garchparam(jj)*Sigma_(:,:,tt-jj);
    end
end
