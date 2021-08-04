function [ nLogL, logLcontr, SigmaE, ScaledScore, varargout ] = ...
    gas_scalar_har_likeRec( param, X, dist )
%
% Michael Stollenwerk
% michael.stollenwerk@live.com
% 11.04.2021

[k,~,T] = size(X);
k_ = k*(k+1)/2;
%% Parameters
intrcpt = ivechchol(param(1:k_));
scoreparam = param(k_ + 1);
dailyparam = param(k_ + 2);
weeklyparam = param(k_ + 3);
monthlyparam = param(k_ + 4);

if strcmp( dist, 'Wish' )
    n = param(k_ + 4 + 1);    
elseif strcmp( dist, 'iWish' )
    n = param(k_ + 4 + 1);  
elseif strcmp( dist, 'tWish' )
    n = param(k_ + 4 + 1); 
    nu = param(k_ + 4 + 2);
elseif strcmp( dist, 'itWish' )
    n = param(k_ + 4 + 1); 
    nu = param(k_ + 4 + 2);
elseif strcmp( dist, 'F' )
    n = param(k_ + 4 + 1); 
    nu = param(k_ + 4 + 2);
elseif strcmp( dist, 'Riesz' )
    n = param(k_ + 4 + 1 : k_ + 4 + k);    
elseif strcmp( dist, 'iRiesz' )
    n = param(k_ + 4 + 1 : k_ + 4 + k);  
elseif strcmp( dist, 'tRiesz' )
    n = param(k_ + 4 + 1 : k_ + 4 + k); 
    nu = param(k_ + 4 + k + 1);
elseif strcmp( dist, 'itRiesz' )
    n = param(k_ + 4 + 1 : k_ + 4 + k); 
    nu = param(k_ + 4 + k + 1);
elseif strcmp( dist, 'FRiesz' )
    n = param(k_ + 4 + 1 : k_ + 4 + k); 
    nu = param(k_ + 4 + k + 1 : k_ + 4 + k + k);        
end

if nargout == 5
    param_out.intrcpt = intrcpt;
    param_out.dailyparam = dailyparam;
    param_out.weeklyparam = weeklyparam;
    param_out.monthlyparam = monthlyparam;
    param_out.n = n;
    
    if strcmp( dist, 'tWish' ) ...
         || strcmp( dist, 'itWish' ) ...
         || strcmp( dist, 'F' ) ...
         || strcmp( dist, 'tRiesz' ) ...
         || strcmp( dist, 'itRiesz' ) ...
         || strcmp( dist, 'FRiesz' )
 
        param_out.nu = nu;
        
    end

    param_out.all = param;
    
    varargout{1} = param_out;
end
%% Data Storage
SigmaE = NaN(k,k,T+1);
ScaledScore = NaN(k,k,T);
logLcontr = NaN(T,1);
%% Likelihood-Recursion
% Initialize recursion at unconditional mean (stationarity assumed).
ini_Sigma = intrcpt./(1 - dailyparam ...
                        - weeklyparam ...
                        - monthlyparam);

for tt=1:T+1
    
    % Recursion                 
    SigmaE(:,:,tt) = intrcpt;
    for jj = 1:22
        if (tt-jj)<1
            if jj == 1
                SigmaE(:,:,tt) = SigmaE(:,:,tt) ...
                                 + dailyparam*ini_Sigma ...
                                 + weeklyparam/5*ini_Sigma ...
                                 + monthlyparam/22*ini_Sigma;
            elseif (1 < jj) && (jj <= 5)
                SigmaE(:,:,tt) = SigmaE(:,:,tt) ...
                                 + weeklyparam/5*ini_Sigma ...
                                 + monthlyparam/22*ini_Sigma;   
            elseif (5 < jj) && (jj <= 22)
                SigmaE(:,:,tt) = SigmaE(:,:,tt) ...
                                 + monthlyparam/22*ini_Sigma;
            end           
        else
            if jj == 1
                SigmaE(:,:,tt) = SigmaE(:,:,tt) ...
                                 + scoreparam*ScaledScore(:,:,tt-jj) ...
                                 + dailyparam*SigmaE(:,:,tt-jj) ...
                                 + weeklyparam/5*SigmaE(:,:,tt-jj) ...
                                 + monthlyparam/22*SigmaE(:,:,tt-jj);
            elseif (1 < jj) && (jj <= 5)
                SigmaE(:,:,tt) = SigmaE(:,:,tt) ...
                                 + weeklyparam/5*SigmaE(:,:,tt-jj) ...
                                 + monthlyparam/22*SigmaE(:,:,tt-jj);   
            elseif (5 < jj) && (jj <= 22)
                SigmaE(:,:,tt) = SigmaE(:,:,tt) ...
                                 + monthlyparam/22*SigmaE(:,:,tt-jj);
            end
        end
    end


    if tt < T+1 
    
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

    end
    
end
%% Log-Likelihood
nLogL = -sum(logLcontr);
