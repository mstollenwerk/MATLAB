function [ nLogL, logLcontr, SigmaE, ScaledScore, varargout ] = ...
    gas_scalar_BEKK_likeRec( param, p, q, X, dist )
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

if strcmp( dist, 'Wish' )
    n = param(k_+ p + q + 1);    
elseif strcmp( dist, 'iWish' )
    n = param(k_+ p + q + 1);  
elseif strcmp( dist, 'tWish' )
    n = param(k_+ p + q + 1); 
    nu = param(k_+ p + q + 2);
elseif strcmp( dist, 'itWish' )
    n = param(k_+ p + q + 1); 
    nu = param(k_+ p + q + 2);
elseif strcmp( dist, 'F' )
    n = param(k_+ p + q + 1); 
    nu = param(k_+ p + q + 2);
elseif strcmp( dist, 'Riesz' )
    n = param(k_+ p + q + 1 : k_+ p + q + k);    
elseif strcmp( dist, 'iRiesz' )
    n = param(k_+ p + q + 1 : k_+ p + q + k);  
elseif strcmp( dist, 'tRiesz' )
    n = param(k_+ p + q + 1 : k_+ p + q + k); 
    nu = param(k_+ p + q + k + 1);
elseif strcmp( dist, 'itRiesz' )
    n = param(k_+ p + q + 1 : k_+ p + q + k); 
    nu = param(k_+ p + q + k + 1);
elseif strcmp( dist, 'FRiesz' )
    n = param(k_+ p + q + 1 : k_+ p + q + k); 
    nu = param(k_+ p + q + k + 1 : k_+ p + q + k + k);        
end

if nargout == 5
    param_out.intrcpt = intrcpt;
    param_out.scoreparam = scoreparam;
    param_out.garch_param = garchparam;
    param_out.n = n;
    
    if exist('nu','var')
        param_out.nu = nu;
    end

    param_out.all = param;
    
    varargout{1} = param_out;
end
%% Data Storage
SigmaE = NaN(k,k,T+1);
ScaledScore = NaN(k,k,T);
logLcontr = NaN(T,1);
%% Recursion
% Initialize recursion at unconditional mean (stationarity assumed).
ini_Sigma = intrcpt./(1 - sum(garchparam));
for tt=1:T
    SigmaE(:,:,tt) = intrcpt;
    for jj = 1:p
        if (tt-jj) <= 0
            SigmaE(:,:,tt) = SigmaE(:,:,tt) + scoreparam(jj)*zeros(k);
        else
            SigmaE(:,:,tt) = SigmaE(:,:,tt) + scoreparam(jj)*ScaledScore(:,:,tt-jj);
        end
    end
    for jj = 1:q
        if (tt-jj) <= 0
            SigmaE(:,:,tt) = SigmaE(:,:,tt) + garchparam(jj)*ini_Sigma;
        else
            SigmaE(:,:,tt) = SigmaE(:,:,tt) + garchparam(jj)*SigmaE(:,:,tt-jj);
        end
    end
    
    try
        % Centering distributions
        if strcmp( dist, 'Wish' )
            Sigma_ = SigmaE(:,:,tt)/n;
        elseif strcmp( dist, 'iWish' )
            Sigma_ = SigmaE(:,:,tt)*(n-k-1);
        elseif strcmp( dist, 'tWish')
            Sigma_ = SigmaE(:,:,tt)*(nu-2)/nu/n;
        elseif strcmp( dist, 'itWish')
            Sigma_ = SigmaE(:,:,tt)*(n-k-1);
        elseif strcmp( dist, 'F' )
            Sigma_ = SigmaE(:,:,tt)*(nu-k-1)/n;
        elseif strcmp( dist, 'Riesz' )
            cholSigE = chol(SigmaE(:,:,tt),'lower');
            Sigma_ = cholSigE/diag(n)*cholSigE';
        elseif strcmp( dist, 'iRiesz' )
            invCholInvSigE = inv(chol(inv(SigmaE(:,:,tt)),'lower'));
            M = matviRieszexpmat(n);
            Sigma_ = invCholInvSigE'/M*invCholInvSigE;
        elseif strcmp( dist, 'tRiesz' )
            cholSigE = chol(SigmaE(:,:,tt),'lower');
            Sigma_ = cholSigE/diag(n)*cholSigE'.*(nu - 2)./nu;
        elseif strcmp( dist, 'itRiesz' )
            invCholInvSigE = inv(chol(inv(SigmaE(:,:,tt)),'lower'));
            M = matviRieszexpmat(n);
            Sigma_ = invCholInvSigE'/M*invCholInvSigE;
        elseif strcmp( dist, 'FRiesz' )
            cholSigE = chol(SigmaE(:,:,tt),'lower');
            M = matvFRieszexpmat(n,nu);
            Sigma_ = cholSigE/M*cholSigE';
        end
        % Likelihood Evaluation    
        if exist('nu','var')
            [logLcontr(tt), score] = logpdf( dist, X(:,:,tt), Sigma_, n, nu );
        else
            [logLcontr(tt), score] = logpdf( dist, X(:,:,tt), Sigma_, n );
        end
        S = ivech(score.Sigma_);
        S = (S + diag(diag(S)))./2; %disregard symmetry in differential
    catch ME
%         ME.message
%         [ME.stack.line]'
%         {ME.stack.name}'
        nLogL = inf;
        logLcontr = NaN;
        SigmaE = NaN;
        param = NaN;        
        return
    end
    
    % Scaled Score 
    if strcmp( dist, 'Wish' )
        S = 2*Sigma_*S*Sigma_;
    elseif strcmp( dist, 'iWish' )
        S = -2/(n-k-1)/n*Sigma_*S*Sigma_;
    elseif strcmp( dist, 'tWish')
        S = 2*Sigma_*S*Sigma_;
    elseif strcmp( dist, 'itWish')
        S = -2/(n-k-1)/n/nu*(nu-2)*Sigma_*S*Sigma_;
    elseif strcmp( dist, 'F' )
        S = 2*Sigma_*S*Sigma_;
    elseif strcmp( dist, 'Riesz' )
        S = 2*Sigma_*S*Sigma_;
    elseif strcmp( dist, 'iRiesz' )
        S = -2*k/sum(diag(M).^(-1))/k*sum(n.^(-1))*Sigma_*S*Sigma_;
    elseif strcmp( dist, 'tRiesz' )
        S = 2*Sigma_*S*Sigma_;
    elseif strcmp( dist, 'itRiesz' )
        S = -2*k/sum(diag(M).^(-1))*sum(n.^(-1))/k/nu*(nu-2)*Sigma_*S*Sigma_;
    elseif strcmp( dist, 'FRiesz' )
        S = 2*Sigma_*S*Sigma_;
    end
    ScaledScore(:,:,tt) = (S+S')./2; % to get rid of small asymmetries.
%         % Scaling as in Opschoor et. al (2018)
%         S = ivech(S);
%             % Opschoor et al. take derivative ignogring that Sigma_ is
%             % symmetric, whereas my scores always take symmetry into
%             % account. To undo that, half the off-diagonal elements of S.
%         S = .5*(ones(p) + eye(p)).*S; 
%         ScaledScore(:,:,tt) =  2/(sum(n.*ones(k,1))/k + 1) * Sigma_(:,:,tt)*S*Sigma_(:,:,tt);   

%         % Scaling as proposed by Creal, Koopman, and Lucas (2013)
%         Score(:,:,tt) = ivech( inv(fisherinfo.Sigma_)*score.Sigma_' );  
end
%% Fcst
for tt=T+1:T+22
    SigmaE(:,:,tt) = intrcpt;
    for jj = 1:p
        if (tt-jj) > T
            SigmaE(:,:,tt) = SigmaE(:,:,tt) + scoreparam(jj)*zeros(k);
        else
            SigmaE(:,:,tt) = SigmaE(:,:,tt) + scoreparam(jj)*ScaledScore(:,:,tt-jj);
        end
    end
    for jj = 1:q
        SigmaE(:,:,tt) = SigmaE(:,:,tt) + garchparam(jj)*SigmaE(:,:,tt-jj);
    end
end
%% Log-Likelihood
nLogL = -sum(logLcontr);
