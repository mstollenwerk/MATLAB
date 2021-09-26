function [eparam, tstats, logL, fit, fcst, optimoutput] = ...
	gas_diag_BEKK_estim_targeting( X, p, q, dist, scalingtype, x0, varargin )
%
%
% Michael Stollenwerk
% michael.stollenwerk@live.com
% 24.09.2021

%% Input Checking
% Will be added later
[k,~,T] = size(X);
k_ = k*(k+1)/2;
%% Optimization
obj_fun = @(param) obj_fun_wrapper(param, X, p, q, dist, scalingtype);
% x0-----------------------------------------------------------------------
if isempty(x0) || (~isempty(x0) && (isinf(obj_fun(x0)) || isnan(obj_fun(x0))))
    if strcmp( dist, 'Wish' )
		x0_df = 2*k;
    elseif strcmp( dist, 'iWish' )
		x0_df = 2*k;  
    elseif strcmp( dist, 'tWish' )
		x0_df = [ 2*k, 5 ];
    elseif strcmp( dist, 'itWish' )
		x0_df = [ 2*k, 5 ];
    elseif strcmp( dist, 'F' )
		x0_df = [ 2*k+3, 2*k+3 ];  % nu > p - 3 for Fisher Info to exist.
    elseif strcmp( dist, 'Riesz' )
		x0_df = ones(1,k).*(2*k);
    elseif strcmp( dist, 'Riesz2' )
		x0_df = ones(1,k).*(2*k);
    elseif strcmp( dist, 'iRiesz' )
		x0_df = ones(1,k).*(2*k);   
    elseif strcmp( dist, 'iRiesz2' )
		x0_df = ones(1,k).*(2*k);   
    elseif strcmp( dist, 'tRiesz' )
		x0_df = [ ones(1,k).*(2*k), 5 ]; 
    elseif strcmp( dist, 'tRiesz2' )
		x0_df = [ ones(1,k).*(2*k), 5 ]; 
    elseif strcmp( dist, 'itRiesz' )
		x0_df = [ ones(1,k).*(2*k), 5 ];
    elseif strcmp( dist, 'itRiesz2' )
		x0_df = [ ones(1,k).*(2*k), 5 ];
    elseif strcmp( dist, 'FRiesz' )
		x0_df = [ ones(1,k).*(2*k+3), ones(1,k).*(2*k+3) ];
    elseif strcmp( dist, 'FRiesz2' )
		x0_df = [ ones(1,k).*(2*k+3), ones(1,k).*(2*k+3) ]; 
    end
    % Candidate Starting Points for Optimization:
    x0_1 = [ones(1,p*k)*.1/p, ones(1,q*k)*.75/q, x0_df]';
    x0_2 = [ones(1,p*k)*.01/p, ones(1,q*k)*.75/q, x0_df]';
    x0_3 = [zeros(1,p*k), ones(1,q*k)*.75/q, x0_df]';
    x0_4 = [zeros(1,p*k+q*k), x0_df]';    
    if ~isinf(obj_fun(x0_1)) && ~isnan(obj_fun(x0_1))
        x0 = x0_1;
    elseif ~isinf(obj_fun(x0_2)) && ~isnan(obj_fun(x0_2))
        x0 = x0_2;
    elseif ~isinf(obj_fun(x0_3)) && ~isnan(obj_fun(x0_3))
        x0 = x0_3;
    elseif ~isinf(obj_fun(x0_4)) && ~isnan(obj_fun(x0_4))
        x0 = x0_4;
    end
end
% Restrictions-------------------------------------------------------------
if strcmp( dist, 'Wish' )
    lb = [-inf(p*k + q*k,1); k-1];
elseif strcmp( dist, 'iWish' )
    lb = [-inf(p*k + q*k,1); k+1];
elseif strcmp( dist, 'F' )
    lb = [-inf(p*k + q*k,1); k-1; k+1];
elseif strcmp( dist, 'tWish' )
    lb = [-inf(p*k + q*k,1); k-1; 2];
elseif strcmp( dist, 'itWish' )
    lb = [-inf(p*k + q*k,1); k+1; 0];    
elseif strcmp( dist, 'Riesz' )
    lb = [-inf(p*k + q*k,1)', 0:k-1]';  
elseif strcmp( dist, 'Riesz2' )
    lb = [-inf(p*k + q*k,1)', fliplr(0:k-1)]'; %see stochstic representation
elseif strcmp( dist, 'iRiesz' )
    lb = [-inf(p*k + q*k,1)', 2:k+1]';
elseif strcmp( dist, 'iRiesz2' )
    lb = [-inf(p*k + q*k,1)', fliplr(2:k+1)]'; %conjectured
elseif strcmp( dist, 'tRiesz' )
    lb = [-inf(p*k + q*k,1)', 0:k-1, 2]'; 
elseif strcmp( dist, 'tRiesz2' )
    lb = [-inf(p*k + q*k,1)', fliplr(0:k-1), 2]'; %see stochstic representation
elseif strcmp( dist, 'itRiesz' )
    lb = [-inf(p*k + q*k,1)', 2:k+1, 2]'; %conjectured
elseif strcmp( dist, 'itRiesz2' )
    lb = [-inf(p*k + q*k,1)', fliplr(2:k+1), 2]'; %conjectured
elseif strcmp( dist, 'FRiesz' )
    lb = [-inf(p*k + q*k,1)', (0:k-1), (k-(1:k)+2)]';
elseif strcmp( dist, 'FRiesz2' )
    lb = [-inf(p*k + q*k,1)', (k-(1:k)), (2:k+1) ]'; %conjectured
end
% A = [ zeros(1,p) ones(1,q) zeros(1, length(x0)-p-q) ];   % Stationarity
% b = 1;                                                   % Stationarity
% Optimization-------------------------------------------------------------
disp('-------------------------------------------------------------------')
disp('-------------------------------------------------------------------')
disp(strcat("Starting gas diag BEKK(",num2str(p),",",num2str(q),")-",dist," estimation, targeting."))
disp('-------------------------------------------------------------------')
disp('-------------------------------------------------------------------')
[eparam,optimoutput] = ...
	my_fmincon_robust(...
		obj_fun, ...
		x0, ...
        lb,[], ...
        [p*k + q*k, numel(x0)]', ...
        1e-2,...
		varargin{:} ...
	);
disp('-------------------------------------------------------------------')
disp('-------------------------------------------------------------------')
disp(strcat("Finished gas diag BEKK(",num2str(p),",",num2str(q),")-",dist," estimation, targeting."))
disp('-------------------------------------------------------------------')
disp('-------------------------------------------------------------------')
% Tstats-------------------------------------------------------------------
%[VCV,A,B,scores,hess,gross_scores] = robustvcv(fun, eparam, 3);
[VCV,scores,gross_scores] = vcv( obj_fun, eparam );
tstats = eparam./sqrt(diag(VCV));

% Output Creation----------------------------------------------------------
[ nLogL, logLcontr, Sigma_, ScaledScore, eparam ] = obj_fun( eparam );

aic = 2*nLogL + 2*(numel(x0)+k_);
bic = 2*nLogL + log(T)*(numel(x0)+k_); % see Yu, Li and Ng (2017) 
logL = struct(...
    'logL', -nLogL,...
    'aic', aic,...
    'bic', bic,...
    'logLcontr', logLcontr...
);

fit = struct( ...
    'Sigma_', Sigma_(:,:,1:T), ...
    'ScaledScore', ScaledScore ...
);
fcst = struct('Sigma_', Sigma_(:,:,T+1:end));

end

function [ nLogL, logLcontr, Sigma_, ScaledScore, param ] = obj_fun_wrapper(param, X, p, q, dist, scalingtype) 
    
    k = size(X,1);
    
    garchparam = reshape(param(p*k+1:p*k+q*k),k,q);
    if any(sum(garchparam,2)>1) % Not entirely sure, based on Engle and Kroner (1995), Section 2.1.
        nLogL = inf;   
        return
    end
    
    meanX = mean(X,3);
    k = size(X,1);
    if strcmp( dist, 'Wish' )
        n = param(p*k + q*k + 1);
        meanSig = meanX/n;
    elseif strcmp( dist, 'iWish' )
        n = param(p*k + q*k + 1);
        meanSig = meanX*(n-k-1); 
    elseif strcmp( dist, 'tWish' )
        n = param(p*k + q*k + 1);
        nu = param(p*k + q*k + 2);
        meanSig = meanX/n/nu*(nu-2);   
    elseif strcmp( dist, 'itWish' )
        n = param(p*k + q*k + 1);
        meanSig = meanX*(n-k-1); 
    elseif strcmp( dist, 'F' )
        n = param(p*k + q*k + 1);
        nu = param(p*k + q*k+2);
        meanSig = meanX/n/nu*(nu-k-1);    
    elseif strcmp( dist, 'Riesz' )
        n = param(p*k + q*k + 1:p*k + q*k+k);
        cholmeanX = chol(meanX,'lower');
        cholmeanSig = cholmeanX/sqrtm(diag(n));
        meanSig = cholmeanSig*cholmeanSig';  
    elseif strcmp( dist, 'Riesz2' )
        n = param(p*k + q*k + 1:p*k + q*k+k);
        cholmeanX = cholU(meanX);
        cholmeanSig = cholmeanX/sqrtm(diag(n));
        meanSig = cholmeanSig*cholmeanSig';        
    elseif strcmp( dist, 'iRiesz' )
        n = param(p*k + q*k + 1:p*k + q*k+k);
        cholmeanX = cholU(meanX);
        cholmeanSig = cholmeanX/sqrtm(matviRieszexpmat(n));
        meanSig = cholmeanSig*cholmeanSig';   
    elseif strcmp( dist, 'iRiesz2' )
        n = param(p*k + q*k + 1:p*k + q*k+k);
        cholmeanX = chol(meanX,'lower');
        cholmeanSig = cholmeanX/sqrtm(matviRiesz2expmat(n));
        meanSig = cholmeanSig*cholmeanSig';   
    elseif strcmp( dist, 'tRiesz' )
        n = param(p*k + q*k + 1:p*k + q*k+k);
        nu = param(p*k + q*k+k + 1);
        cholmeanX = chol(meanX,'lower');
        cholmeanSig = cholmeanX/sqrtm(diag(n)*nu/(nu-2));
        meanSig = cholmeanSig*cholmeanSig';
    elseif strcmp( dist, 'tRiesz2' )
        n = param(p*k + q*k + 1:p*k + q*k+k);
        nu = param(p*k + q*k+k + 1);
        cholmeanX = cholU(meanX);
        cholmeanSig = cholmeanX/sqrtm(diag(n)*nu/(nu-2));
        meanSig = cholmeanSig*cholmeanSig';   
    elseif strcmp( dist, 'itRiesz' )
        n = param(p*k + q*k + 1:p*k + q*k+k);
        cholmeanX = cholU(meanX);
        cholmeanSig = cholmeanX/sqrtm(matviRieszexpmat(n));
        meanSig = cholmeanSig*cholmeanSig';   
    elseif strcmp( dist, 'itRiesz2' )
        n = param(p*k + q*k + 1:p*k + q*k+k);
        cholmeanX = chol(meanX,'lower');
        cholmeanSig = cholmeanX/sqrtm(matviRiesz2expmat(n));
        meanSig = cholmeanSig*cholmeanSig';
    elseif strcmp( dist, 'FRiesz' )
        n = param(p*k + q*k + 1:p*k + q*k+k);
        nu = param(p*k + q*k+k + 1:p*k + q*k+k+k);
        cholmeanX = chol(meanX,'lower');
        cholmeanSig = cholmeanX/sqrtm(matvFRieszexpmat(n,nu));
        meanSig = cholmeanSig*cholmeanSig';
    elseif strcmp( dist, 'FRiesz2' )
        n = param(p*k + q*k + 1:p*k + q*k+k);
        nu = param(p*k + q*k+k + 1:p*k + q*k+k+k);
        cholmeanX = cholU(meanX);
        cholmeanSig = cholmeanX/sqrtm(matvFRiesz2expmat(n,nu));
        meanSig = cholmeanSig*cholmeanSig';
    end
    
    Intrcpt = meanSig;
    for ii = 1:q
        B = diag(sqrt(garchparam(:,ii)));
        Intrcpt = Intrcpt - B*meanSig*B';
    end
    try
        vechcholIntrcpt = vechchol(Intrcpt);
    catch
        nLogL = inf;   
        return
    end
    
    [ nLogL, logLcontr, Sigma_, ScaledScore, param ] = gas_diag_BEKK_likeRec( ...
        [vechcholIntrcpt; param], ...
        p, ...
        q, ...
        X, ...
        dist, ...
        scalingtype ...
    );

end