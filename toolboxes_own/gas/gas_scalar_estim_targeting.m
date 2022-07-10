function [eparam, tstats, logL, fit, fcst, optimoutput] = ...
	gas_scalar_estim_targeting( R, dist, x0, varargin )
%
%
% Michael Stollenwerk
% michael.stollenwerk@live.com
% 30.06.2022

%% Input Checking
% Will be added later
[p,~,T] = size(R);
p_ = p*(p+1)/2;
%% Optimization
obj_fun = @(param, perm_) obj_fun_wrapper(param, R(perm_,perm_,:), dist);
% x0-----------------------------------------------------------------------
if ~isempty(x0) && (isinf(obj_fun(x0,1:p)) || isnan(obj_fun(x0,1:p)))
    error('User supplied x0 invalid.')
end
% The lower bounds below are for existence of the distributions. For
% existence of the expected value matrix Sigma_, restrictions on the dfs
% are stricter. The inverse Wishart, for example, exists for nu > p-1, but
% its mean only exists for nu > p+1.
lb_df_barlett = [0, -inf, 1e-5];
lb_df_barlettL = [zeros(1,p), -inf, 1e-5];
lb_df_barlettU = [zeros(1,p), -inf, 1e-5];
lb_df_chi = [0, -inf, 1e-5];
ub_df_barlett = [inf(1,2), .9999];
ub_df_barlettL = [inf(1,p+1), .999];
ub_df_barlettU = [inf(1,p+1), .999];
ub_df_chi = [inf(1,2), .9999];
x0_df_barlett = [2*p.*.01, 0.00001, .99];
x0_df_barlettL = [2*p*ones(1,p).*.01, 0.00001, .99];
x0_df_barlettU = [2*p*ones(1,p).*.01, 0.00001, .99];
x0_df_chi = [5, 0.00001, .99];

if strcmp( dist, 'Wish' )
    x0_df = x0_df_barlett;
    lb_df = lb_df_barlett;
    ub_df = ub_df_barlett;
%     Aeq = [ zeros(4,6) eye(4) ];
%     beq = zeros(4,1);       
elseif strcmp( dist, 'iWish' )
    x0_df = x0_df_barlett;  
    lb_df = lb_df_barlett;
    ub_df = ub_df_barlett;
%     Aeq = [ zeros(4,5*p+1) eye(4) ];
%     beq = zeros(4,1);    
elseif strcmp( dist, 'tWish' )
    x0_df = [ x0_df_barlett, x0_df_chi ];
    lb_df = [ lb_df_barlett, lb_df_chi ];
    ub_df = [ ub_df_barlett, ub_df_chi ];
%     Aeq = [ zeros(4,6) eye(4) zeros(4,5) ;
%             zeros(4,11) eye(4)];
%     beq = zeros(8,1);    
elseif strcmp( dist, 'itWish' )
    x0_df = [ x0_df_chi, x0_df_barlett ];
    lb_df = [ lb_df_chi, lb_df_barlett ];
    ub_df = [ ub_df_chi, ub_df_barlett ];
%     Aeq = [ zeros(4,6) eye(4) zeros(4,5) ;
%             zeros(4,11) eye(4)];
%     beq = zeros(8,1);      
elseif strcmp( dist, 'F' )
    x0_df = [ x0_df_barlett, x0_df_barlett ];   
    lb_df = [ lb_df_barlett, lb_df_barlett ];
    ub_df = [ ub_df_barlett, ub_df_barlett ];
%     Aeq = [ zeros(4,6) eye(4) zeros(4,5) ;
%             zeros(4,11) eye(4)];
%     beq = zeros(8,1);        
elseif strcmp( dist, 'Riesz' )
    x0_df = x0_df_barlettL;
    lb_df = lb_df_barlettL;
    ub_df = ub_df_barlettL;
elseif strcmp( dist, 'iRiesz2' )
    x0_df = x0_df_barlettU;
    lb_df = lb_df_barlettU;
    ub_df = ub_df_barlettU;
elseif strcmp( dist, 'tRiesz' )
    x0_df = [ x0_df_barlettL, x0_df_chi ]; 
    lb_df = [ lb_df_barlettL, lb_df_chi ];
    ub_df = [ ub_df_barlettL, ub_df_chi ];
elseif strcmp( dist, 'itRiesz2' )
    x0_df = [ x0_df_chi, x0_df_barlettU ];
    lb_df = [ lb_df_chi, lb_df_barlettU];
    ub_df = [ ub_df_chi, ub_df_barlettU];
elseif strcmp( dist, 'FRiesz' )
    x0_df = [ x0_df_barlettL, x0_df_barlettU ];
    lb_df = [ lb_df_barlettL, lb_df_barlettU];
    ub_df = [ ub_df_barlettL, ub_df_barlettU];
elseif strcmp( dist, 'iFRiesz2' )
    x0_df = [ x0_df_barlettL, x0_df_barlettU ];
    lb_df = [ lb_df_barlettL, lb_df_barlettU ];
    ub_df = [ ub_df_barlettL, ub_df_barlettU ];
end
if isempty(x0)
    x0 = [0.0005, 0.0005, 0.7, 0.1, 0.05, x0_df];
end

% A = [ zeros(1,p) ones(1,q) zeros(1, length(x0)-p-q) ];   % Stationarity
% b = 1;                                                   % Stationarity
lb = [zeros(1,5), lb_df];
ub = [inf, inf, ones(1,3), ub_df];

if contains(dist,'Riesz')
    perm_optim = 1;
else
    perm_optim = 0;
end
if ~isempty(varargin) && strcmp(varargin{1},'noperm')
    perm_optim = 0;
    varargin = varargin(2:end);
end
if perm_optim
    [eparam,optimoutput] = ...
        fmincon_Rieszperm(...
            p, ...
            obj_fun, ...
            x0', ...
            [],[],[],[], ...
            lb,ub,[], ...
            varargin{:} ...
        );
    obj_fun_opt_perm = @(param) obj_fun(param,optimoutput.perm_);
else
    [eparam,optimoutput] = ...
        my_fmincon(...
            @(param) obj_fun(param,1:p), ...
            x0', ...
            [],[],[],[], ...
            lb,ub,[], ...
            varargin{:} ...
        );
    obj_fun_opt_perm = @(param) obj_fun(param,1:p);
    optimoutput.perm_ = 1:p;
end
    
% Output Creation----------------------------------------------------------
[ nLogL, logLcontr, dyn, S, eparam ] = obj_fun_opt_perm( eparam );
eparam.perm_ = optimoutput.perm_;

obj_fun_no_targeting = @(param) gas_scalar_likeRec( ...
        [eparam.all(1:p_); param], ...
        R(eparam.perm_,eparam.perm_,:), ...
        dist ...
    );
%[VCV,A,B,scores,hess,gross_scores] = robustvcv(fun, eparam, 3);
[VCV,scores,gross_scores] = vcv( obj_fun_no_targeting, eparam.all(p_+1:end) );
tstats = eparam.all(p_+1:end)./sqrt(diag(VCV));


aic = 2*nLogL + 2*(numel(x0)+p_);
bic = 2*nLogL + log(T)*(numel(x0)+p_); % see Yu, Li and Ng (2017) 
logL_detR_part = -(p+1)/2*sum(logdet3d(R));
logL = struct(...
    'logL', -nLogL,...
    'logL_detR_part', logL_detR_part, ...
    'aic', aic,...
    'bic', bic,...
    'logLcontr', logLcontr...
);

fit = struct( ...
    'Sig', dyn.Sig(:,:,1:T), ...
    'ScaledScoreSig', S.Sig(:,:,1:T) ...
);
fcst = struct( ...
    'Sig', dyn.Sig(:,:,T+1:end) ...
);
if isfield(S,'n')
    fit.ScaledScoreN = S.n(:,1:T);
    fit.n = dyn.n(:,1:T);
    fcst.n = dyn.n(:,T+1:end);
end
if isfield(S,'nu')
    fit.ScaledScoreNu = S.nu(:,1:T);
    fit.nu = dyn.nu(:,1:T);
    fcst.nu = dyn.nu(:,T+1:end);
end

end

function [ nLogL, logLcontr, dyn, S, param_out ] = ...
    obj_fun_wrapper(param, R, dist) 
    
    p = size(R,1);
    
    if size(param,1)<size(param,2)
        param = param';
    end
    
    persistenceSig = param(3) + param(4) + param(5);
    if persistenceSig >= 1
        nLogL = inf;
        return
    end
    
    meanSig = mean(R,3);
    
    vechcholIntrcptSig = vechchol( (1-persistenceSig)*meanSig );  
    
    [ nLogL, logLcontr, dyn, S, param_out ] = gas_scalar_likeRec( ...
        [vechcholIntrcptSig; param], ...
        R, ...
        dist ...
    );

end