function [eparam, tstats, logL, fit, fcst, weights, optimoutput] = ...
	caicov_scalar_BEKK_generalized_estim_targeting( X, p, q, x0, varargin )
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
obj_fun = @(param) obj_fun_wrapper(param, X, p, q);
% x0-----------------------------------------------------------------------
if ~isempty(x0) && (isinf(obj_fun(x0)) || isnan(obj_fun(x0)))
    error('User supplied x0 invalid.')
end
if isempty(x0)
%     elseif strcmp( dist, 'tRiesz' )
% 		x0_df = [ ones(1,k).*(2*k), 5 ]; 
    x0_df = [ ones(1,k).*(2*k), 5, ones(1,k).*(2*k+3), ones(1,k).*(2*k+3) ];
    % Candidate Starting Points for Optimization:
%     x0 = [ones(1,p)*.1/p, ones(1,q)*.75/q, 0, 0, 0, 0, 0, x0_df;
%           ones(1,p)*.01/p, ones(1,q)*.75/q, .5, .5, 0, .9, .01, x0_df;
%           zeros(1,p), ones(1,q)*.75/q, .5, .5, 0, .9, .01, x0_df;
%           zeros(1,p+q), .5, .5, 0, .9, .01, x0_df]';  
    x0 = [ones(1,p)*.1/p, ones(1,q)*.75/q, 0, 0, 0, 0, x0_df;
          ones(1,p)*.01/p, ones(1,q)*.75/q, .5, .5, 0, .01, x0_df;
          zeros(1,p), ones(1,q)*.75/q, .5, .5, 0, .01, x0_df;
          zeros(1,p+q), .5, .5, 0, .9, x0_df]';    
end
% Restrictions-------------------------------------------------------------
%     lb = [-inf(p+q+3,1)', 0:k-1, 2]'; 
% lb = [0, 0, -1, -inf, -inf, -inf, -inf, 2:k+1, 2, (0:k-1), (k-(1:k)+2)]'; %conjectured
% ub = [1, 1, 1, inf, inf, inf, inf, inf(1,k+1+k+k)]';
lb = [0, 0, 0, 0, 0, -1, 2:k+1, 2, (0:k-1), (k-(1:k)+2)]'; %conjectured
ub = [1, 1, 1, 1, 1, 1, inf(1,k+1+k+k)]';
% Optimization-------------------------------------------------------------
for ii = 1:size(x0,2) % Looping through candidate starting points, until one works.
    try
        disp('-------------------------------------------------------------------')
        disp('-------------------------------------------------------------------')
        disp(strcat("Optimization Starting Point for Recursion Parameters: p = (", num2str(x0(1:p,ii)'), "), q = (", num2str(x0(p+1:p+q,ii)'), ")."))
        disp(strcat("Starting caicov scalar BEKK generalized(",num2str(p),",",num2str(q),")- estimation, targeting."))
        disp('-------------------------------------------------------------------')
        disp('-------------------------------------------------------------------')
        [eparam,optimoutput] = ...
            my_fmincon_robust(...
                obj_fun, ...
                x0(:,ii), ...
                lb,ub, ...
                [p+q, p+q+4, p+q+4+k+k+k+1]', ...
                1e-2,...
                varargin{:} ...
            );
%         [eparam,optimoutput] = ...
%             my_fmincon(...
%                 obj_fun, ...
%                 x0(:,ii), ...
%                 [],[],[],[], ...
%                 lb,ub, ...
%                 [], ...
%                 varargin{:} ...
%             );        
        disp('-------------------------------------------------------------------')
        disp('-------------------------------------------------------------------')
        disp(strcat("Finished caicov scalar BEKK generalized(",num2str(p),",",num2str(q),")-estimation, targeting."))
        disp('-------------------------------------------------------------------')
        disp('-------------------------------------------------------------------')
        break
    catch ME
        ME.message
        [ME.stack.line]'
        {ME.stack.name}'
        disp('-------------------------------------------------------------------')
        disp('-------------------------------------------------------------------')
        disp(strcat("Optimization Starting Point for Recursion Parameters: p = (", num2str(x0(1:p,ii)'), "), q = (", num2str(x0(p+1:p+q,ii)'), ") didn't work."))
        disp('-------------------------------------------------------------------')
        disp('-------------------------------------------------------------------')        
    end
end        
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
[ nLogL, logLcontr, Sigma_, weights, eparam ] = obj_fun( eparam );

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

function [ nLogL, logLcontr, Sigma_, weights, param, fitplot ] = obj_fun_wrapper(param, X, p, q) 

    if sum(param(1:p+q)) >= 1
        nLogL = inf;
        return
    end
    
    vechcholSig = vechchol( mean(X,3)*(1-sum(param(1:p+q))) );
    
    [ nLogL, logLcontr, Sigma_, weights, param, fitplot ] = caicov_scalar_BEKK_generalized_likeRec( ...
        [vechcholSig; param], ...
        p, ...
        q, ...
        X ...
    );

end