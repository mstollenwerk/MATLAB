function [ R, Sigma_, S, logLcontr, varargout ] = gas_scalar_BEKK11_rnd( param, dist, k, T )
%
% Michael Stollenwerk
% michael.stollenwerk@live.com
% 22.11.2022

k_ = k*(k+1)/2;
%% Parameters
intrcpt = ivechchol(param(1:k_));
scoreparam = param(k_ + 1);
garchparam = param(k_ + 2);
dfs = param(k_+3:end);
%% Data Storage
R = NaN(k,k,T+1e4);
Sigma_ = NaN(k,k,T+1e4);
S = NaN(k,k,T+1e4);
logLcontr = NaN(T+1e4,1);
%% Recursion
for tt=1:T+1e4
    
    if tt == 1
        Sigma_(:,:,tt) = intrcpt./(1 - sum(garchparam));
    else
        Sigma_(:,:,tt) = intrcpt + scoreparam*S(:,:,tt-1) ...
                                 + garchparam*Sigma_(:,:,tt-1);
    end
    
    R(:,:,tt) = matvsrnd( dist, Sigma_(:,:,tt), dfs, 1 );

    % Likelihood Evaluation         
    [~, logLcontr(tt), score, ~, param_dist] = ...
        matvsLogLike( dist, Sigma_(:,:,tt), dfs, R(:,:,tt) );   
    
    % Scaled Score 
    if any(strcmp( dist, {'F','iFRiesz2','FRiesz'} ))
        n = mean(dfs(1:length(dfs)/2));
        nu = mean(dfs(length(dfs)/2+1:end));
        c_1 = (n^2*(nu-k-2) + 2*n)/((nu-k)*(nu-k-1)*(nu-k-3));
        c_2 = (n*(nu-k-2)+n^2+n)/((nu-k)*(nu-k-1)*(nu-k-3));
        c_4 = (n-k-1)/((n+nu-1)*(n+nu+2))*((n-k-2+1/(nu+n))*c_2-(1+(n-k-1)/(n+nu))*c_1);
        c_3 = (n-k-1)/(n+nu)*((n-k-2)*c_2 - c_1)-(n+nu+1)*c_4;
        alpha = (nu-(n+nu)*(c_3+c_4))/2;
        c = (n+nu)*c_4/2;
    elseif any(strcmp( dist, {'Wish','Riesz'} ))
        n = mean(dfs);
        alpha = n/2;
        c = 0;
    elseif any(strcmp( dist, {'iWish','iRiesz2'} ))
        nu = mean(dfs);
        alpha = nu/2;
        c = 0;
    elseif any(strcmp( dist, {'tWish','tRiesz'} ))
        n = mean(dfs(1:end-1));
        nu = dfs(end);
        alpha = n/2*(nu + k*n)/(nu + k*n + 2);
        c = n/2*n/(nu + k*n + 2);
    elseif any(strcmp( dist, {'itWish','itRiesz2'} ))
        n = dfs(1);
        nu = mean(dfs(2:end));
        alpha = nu/2*(n + k*nu)/(n + k*nu + 2);
        c = nu/2*nu/(n + k*nu + 2);
    end
    beta = alpha*c/(alpha + c*k);    
        
    SigDel = Sigma_(:,:,tt)*score.SigmaNonSym;
    SigDelSig = SigDel*Sigma_(:,:,tt);
    S(:,:,tt) = alpha/2*SigDelSig+SigDelSig' ...
                + beta*trace(SigDel)*Sigma_(:,:,tt);
    
end
R = R(:,:,1e4+1:end);
Sigma_ = Sigma_(:,:,1e4+1:end);
S = S(:,:,1e4+1:end);
logLcontr = logLcontr(1e4+1:end);
%% Parameter Output
if nargout >= 5
    param_out.intrcpt = intrcpt;
    param_out.scoreparam = scoreparam;
    param_out.garch_param = garchparam;
    
    if isfield(param_dist,'n')
        param_out.n = param_dist.n;
    end
    
    if isfield(param_dist,'nu')
        param_out.nu = param_dist.nu;
    end

    param_out.all = param;
    
    varargout{1} = param_out;
end
