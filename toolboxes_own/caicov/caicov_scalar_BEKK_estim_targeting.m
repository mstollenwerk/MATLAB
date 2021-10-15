function [eparam, tstats, logL, fit, fcst, optimoutput] = ...
	caicov_scalar_BEKK_estim_targeting( X, p, q, dist, x0, varargin )
%
%
% Michael Stollenwerk
% michael.stollenwerk@live.com
% 05.04.2021

%% Input Checking
% Will be added later
[k,~,T] = size(X);
k_ = k*(k+1)/2;

%% Optimization
obj_fun = @(param) obj_fun_wrapper(param, X, p, q, dist);
% x0-----------------------------------------------------------------------
if ~isempty(x0) && (isinf(obj_fun(x0)) || isnan(obj_fun(x0)))
    error('User supplied x0 invalid.')
end
if isempty(x0)
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
    x0_1 = [ones(1,p)*.1/p, ones(1,q)*.75/q, x0_df]';
    x0_2 = [ones(1,p)*.01/p, ones(1,q)*.75/q, x0_df]';
    x0_3 = [zeros(1,p), ones(1,q)*.75/q, x0_df]';
    x0_4 = [zeros(1,p+q), x0_df]';    
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
    lb = [-inf(p+q,1); k-1];
elseif strcmp( dist, 'iWish' )
    lb = [-inf(p+q,1); k+1];
elseif strcmp( dist, 'F' )
    lb = [-inf(p+q,1); k-1; k+1];
elseif strcmp( dist, 'tWish' )
    lb = [-inf(p+q,1); k-1; 2];
elseif strcmp( dist, 'itWish' )
    lb = [-inf(p+q,1); k+1; 0];    
elseif strcmp( dist, 'Riesz' )
    lb = [-inf(p+q,1)', 0:k-1]';  
elseif strcmp( dist, 'Riesz2' )
    lb = [-inf(p+q,1)', fliplr(0:k-1)]'; %see stochstic representation
elseif strcmp( dist, 'iRiesz' )
    lb = [-inf(p+q,1)', 2:k+1]';
elseif strcmp( dist, 'iRiesz2' )
    lb = [-inf(p+q,1)', fliplr(2:k+1)]'; %conjectured
elseif strcmp( dist, 'tRiesz' )
    lb = [-inf(p+q,1)', 0:k-1, 2]'; 
elseif strcmp( dist, 'tRiesz2' )
    lb = [-inf(p+q,1)', fliplr(0:k-1), 2]'; %see stochstic representation
elseif strcmp( dist, 'itRiesz' )
    lb = [-inf(p+q,1)', 2:k+1, 2]'; %conjectured
elseif strcmp( dist, 'itRiesz2' )
    lb = [-inf(p+q,1)', fliplr(2:k+1), 2]'; %conjectured
elseif strcmp( dist, 'FRiesz' )
    lb = [-inf(p+q,1)', (0:k-1), (k-(1:k)+2)]';
elseif strcmp( dist, 'FRiesz2' )
    lb = [-inf(p+q,1)', (k-(1:k)+2), (0:k-1) ]'; %conjectured
end
% A = [ ones(1,p+q) zeros(1, length(x0)-p-q) ];   % Stationarity
% b = 1;                                          % Stationarity
% Optimization-------------------------------------------------------------
disp('-------------------------------------------------------------------')
disp('-------------------------------------------------------------------')
disp(strcat("Optimization Starting Point for Recursion Parameters: p = (", num2str(x0(1:p)'), "), q = (", num2str(x0(p+1:p+q)'), ")."))
disp(strcat("Starting caicov scalar BEKK(",num2str(p),",",num2str(q),")-",dist," estimation, targeting."))
disp('-------------------------------------------------------------------')
disp('-------------------------------------------------------------------')
[eparam,optimoutput] = ...
	my_fmincon_robust(...
		obj_fun, ...
		x0, ...
        lb,[], ...
        [p+q, numel(x0)]', ...
        1e-2,...
		varargin{:} ...
 	);
disp('-------------------------------------------------------------------')
disp('-------------------------------------------------------------------')
disp(strcat("Finished caicov scalar BEKK(",num2str(p),",",num2str(q),")-",dist," estimation, targeting."))
disp('-------------------------------------------------------------------')
disp('-------------------------------------------------------------------')
% [eparam,optimoutput] = ...
% 	my_fmincon(...
% 		obj_fun, ...
% 		x0, ...
%         A,b,[],[],lb,[],[], ...
% 		varargin{:} ...
% 	);
% while strcmp(optimoutput.message(50:111),'fmincon stopped because the gradient calculation is undefined.')
%     x0_new = eparam;
%     x0_new(1:q) = x0_new(1:q)*.9;
%     [eparam,optimoutput] = ...
%         my_fmincon(...
%             obj_fun, ...
%             x0_new, ...
%             A,b,[],[],lb,[],[], ...
%             varargin{:} ...
%         );
% end
% Tstats-------------------------------------------------------------------
%[VCV,A,B,scores,hess,gross_scores] = robustvcv(fun, eparam, 3);
[VCV,scores,gross_scores] = vcv( obj_fun, eparam );
tstats = eparam./sqrt(diag(VCV));

% Output Creation----------------------------------------------------------
[ nLogL, logLcontr, Sigma_, eparam ] = obj_fun( eparam );

aic = 2*nLogL + 2*(numel(x0)+k_);
bic = 2*nLogL + log(T)*(numel(x0)+k_); % see Yu, Li and Ng (2017) 
logL = struct(...
    'logL', -nLogL,...
    'aic', aic,...
    'bic', bic,...
    'logLcontr', logLcontr...
);

fit = struct( ...
    'Sigma_', Sigma_(:,:,1:T) ...
);
fcst = struct('Sigma_', Sigma_(:,:,T+1:end));

end

function [ nLogL, logLcontr, Sigma_, param, fitplot ] = obj_fun_wrapper(param, X, p, q, dist) 

    if sum(param(1:p+q)) >= 1
        nLogL = inf;
        return
    end

    vechcholSig = vechchol( mean(X,3)*(1-sum(param(1:p+q))) );
    
    [ nLogL, logLcontr, Sigma_, param, fitplot ] = caicov_scalar_BEKK_likeRec( ...
        [vechcholSig; param], ...
        p, ...
        q, ...
        X, ...
        dist ...
    );

end