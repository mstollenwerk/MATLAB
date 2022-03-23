function [eparam, tstats, logL, fit, fcst, optimoutput] = ...
	gas_scalar_BEKK_estim_targeting( R, p, q, dist, scalingtype, x0, varargin )
%
%
% Michael Stollenwerk
% michael.stollenwerk@live.com
% 25.03.2021

%% Input Checking
% Will be added later
[k,~,T] = size(R);
k_ = k*(k+1)/2;
%% Optimization
obj_fun = @(param, perm_) obj_fun_wrapper(param, R(perm_,perm_,:), p, q, dist, scalingtype);
% x0-----------------------------------------------------------------------
if ~isempty(x0) && (isinf(obj_fun(x0,1:k)) || isnan(obj_fun(x0,1:k)))
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
%     elseif strcmp( dist, 'Riesz2' )
% 		x0_df = ones(1,k).*(2*k);
%     elseif strcmp( dist, 'iRiesz' )
% 		x0_df = ones(1,k).*(2*k);   
    elseif strcmp( dist, 'iRiesz2' )
		x0_df = ones(1,k).*(2*k);   
    elseif strcmp( dist, 'tRiesz' )
		x0_df = [ ones(1,k).*(2*k), 5 ]; 
%     elseif strcmp( dist, 'tRiesz2' )
% 		x0_df = [ ones(1,k).*(2*k), 5 ]; 
%     elseif strcmp( dist, 'itRiesz' )
% 		x0_df = [ ones(1,k).*(2*k), 5 ];
    elseif strcmp( dist, 'itRiesz2' )
		x0_df = [ ones(1,k).*(2*k), 5 ];
    elseif strcmp( dist, 'FRiesz' )
		x0_df = [ ones(1,k).*(2*k+3), ones(1,k).*(2*k+3) ];
%     elseif strcmp( dist, 'FRiesz2' )
% 		x0_df = [ ones(1,k).*(2*k+3), ones(1,k).*(2*k+3) ]; 
    end
    % Candidate Starting Points for Optimization:
    x0 = [zeros(1,p+q), x0_df]'; % This x0 should always work, so you can get rid of all the other ones and the try/catch around the optimization.
%           ones(1,p)*.05/p, ones(1,q)*.98/q, x0_df]';    
end
% Optimization Specs ------------------------------------------------------
if strcmp( dist, 'Wish' )
    lb = [-0.001; -inf(p+q-1,1); k-1];
    perm_optim = 0;
elseif strcmp( dist, 'iWish' )
    lb = [-0.001; -inf(p+q-1,1); k+1];
    perm_optim = 0;
elseif strcmp( dist, 'F' )
    lb = [-0.001; -inf(p+q-1,1); k-1; k+1];
    perm_optim = 0;
elseif strcmp( dist, 'tWish' )
    lb = [-0.001; -inf(p+q-1,1); k-1; 2];
    perm_optim = 0;
elseif strcmp( dist, 'itWish' )
    lb = [-0.001; -inf(p+q-1,1); k+1; 0];    
    perm_optim = 0;
elseif strcmp( dist, 'Riesz' )
    lb = [-0.001, -inf(p+q-1,1)', 0:k-1]';  
    perm_optim = 1;
% elseif strcmp( dist, 'Riesz2' )
%     lb = [-0.001, -inf(p+q-1,1)', fliplr(0:k-1)]'; %see stochstic representation
% elseif strcmp( dist, 'iRiesz' )
%     lb = [-0.001, -inf(p+q-1,1)', 2:k+1]';
elseif strcmp( dist, 'iRiesz2' )
    lb = [-0.001, -inf(p+q-1,1)', fliplr(2:k+1)]'; %conjectured
    perm_optim = 1;
elseif strcmp( dist, 'tRiesz' )
    lb = [-0.001, -inf(p+q-1,1)', 0:k-1, 2]'; 
    perm_optim = 1;
% elseif strcmp( dist, 'tRiesz2' )
%     lb = [-0.001, -inf(p+q-1,1)', fliplr(0:k-1), 2]'; %see stochstic representation
% elseif strcmp( dist, 'itRiesz' )
%     lb = [-0.001, -inf(p+q-1,1)', 2:k+1, 2]'; %conjectured
elseif strcmp( dist, 'itRiesz2' )
    lb = [-0.001, -inf(p+q-1,1)', fliplr(2:k+1), 2]'; %conjectured
    perm_optim = 1;
elseif strcmp( dist, 'FRiesz' )
    lb = [-0.001, -inf(p+q-1,1)', (0:k-1), (k-(1:k)+2)]'; %Note that scoreparam(1) has a high LB. Unfortunatelx FRiesz estimation seems very hard.
    perm_optim = 1;
% elseif strcmp( dist, 'FRiesz2' )
%     lb = [-0.001, -inf(p+q-1,1)', (k-(1:k)), (2:k+1) ]'; %conjectured %Note that scoreparam(1) has a high LB. Unfortunatelx FRiesz estimation seems very hard.
end
% A = [ zeros(1,p) ones(1,q) zeros(1, length(x0)-p-q) ];   % Stationarity
% b = 1;                                                   % Stationarity
% Optimization-------------------------------------------------------------
for ii = 1:size(x0,2) % Looping through candidate starting points, until one works.
    try
        %disp('-------------------------------------------------------------------')
        %disp('-------------------------------------------------------------------')
        %disp(strcat("Optimization Starting Point for Recursion Parameters: p = (", num2str(x0(1:p,ii)'), "), q = (", num2str(x0(p+1:p+q,ii)'), ")."))
        %disp(strcat("Starting gas scalar BEKK(",num2str(p),",",num2str(q),")-",dist," estimation, targeting."))
        %disp('-------------------------------------------------------------------')
        %disp('-------------------------------------------------------------------')
        if perm_optim
            [eparam,optimoutput] = ...
                fmincon_Rieszperm(...
                    k, ...
                    obj_fun, ...
                    x0(:,ii), ...
                    [],[],[],[], ...
                    lb,[],[], ...
                    varargin{:} ...
                );
            obj_fun_opt_perm = @(param) obj_fun(param,optimoutput.perm_);
        else
            [eparam,optimoutput] = ...
                my_fmincon(...
                    @(param) obj_fun(param,1:k), ...
                    x0(:,ii), ...
                    [],[],[],[], ...
                    lb,[],[], ...
                    varargin{:} ...
                );
            obj_fun_opt_perm = @(param) obj_fun(param,1:k);
            optimoutput.perm_ = 1:k;
        end
        %disp('-------------------------------------------------------------------')
        %disp('-------------------------------------------------------------------')
        %disp(strcat("Finished gas scalar BEKK(",num2str(p),",",num2str(q),")-",dist," estimation, targeting."))
        %disp('-------------------------------------------------------------------')
        %disp('-------------------------------------------------------------------')
        break
    catch ME
        ME.message
        [ME.stack.line]'
        {ME.stack.name}'        
        %disp('-------------------------------------------------------------------')
        %disp('-------------------------------------------------------------------')
        %disp(strcat("Optimization Starting Point for Recursion Parameters: p = (", num2str(x0(1:p,ii)'), "), q = (", num2str(x0(p+1:p+q,ii)'), ") didn't work."))
        %disp('-------------------------------------------------------------------')
        %disp('-------------------------------------------------------------------')        
    end
end
    
% Output Creation----------------------------------------------------------
[ nLogL, logLcontr, SigmaE, ScaledScore, eparam ] = obj_fun_opt_perm( eparam );
eparam.perm_ = optimoutput.perm_;

obj_fun_no_targeting = @(param) gas_scalar_BEKK_likeRec( ...
        [eparam.all(1:k_); param], ...
        p, ...
        q, ...
        R(eparam.perm_,eparam.perm_,:), ...
        dist, ...
        scalingtype ...
    );
%[VCV,A,B,scores,hess,gross_scores] = robustvcv(fun, eparam, 3);
[VCV,scores,gross_scores] = vcv( obj_fun_no_targeting, eparam.all(k_+1:end) );
tstats = eparam.all(k_+1:end)./sqrt(diag(VCV));


aic = 2*nLogL + 2*(numel(x0)+k_);
bic = 2*nLogL + log(T)*(numel(x0)+k_); % see Yu, Li and Ng (2017) 
logL_detR_part = -(k+1)/2*sum(logdet3d(R));
logL = struct(...
    'logL', -nLogL,...
    'logL_detR_part', logL_detR_part, ...
    'aic', aic,...
    'bic', bic,...
    'logLcontr', logLcontr...
);

fit = struct( ...
    'SigmaE', SigmaE(:,:,1:T), ...
    'ScaledScore', ScaledScore ...
);
fcst = struct('SigmaE', SigmaE(:,:,T+1:end));

end

function [ nLogL, logLcontr, SigmaE, ScaledScore, param_out, fitplot ] = obj_fun_wrapper(param, R, p, q, dist, scalingtype) 

    if sum(param(p+1:p+q)) >= 1
        nLogL = inf;   
        return
    end
    
    meanSig = mean(R,3);
    
    vechcholIntrcpt = vechchol( meanSig*(1-sum(param(p+1:p+q))) );  
    
    [ nLogL, logLcontr, SigmaE, ScaledScore, param_out, fitplot ] = gas_scalar_BEKK_likeRec( ...
        [vechcholIntrcpt; param], ...
        p, ...
        q, ...
        R, ...
        dist, ...
        scalingtype ...
    );

end