function [eparam, tstats, logL, fit, fcst, optimoutput] = ...
	matvMF2_scalar_BEKK_estim( R, dist, x0, varargin )
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
obj_fun = @(param, perm_) matvMF2_scalar_BEKK_likeRec(param, R(perm_,perm_,:), dist);
% x0-----------------------------------------------------------------------
if ~isempty(x0) && (isinf(obj_fun(x0,1:p)) || isnan(obj_fun(x0,1:p)))
    error('User supplied x0 invalid.')
end
% The lower bounds below are for existence of the distributions. For
% existence of the expected value matrix Sigma_, restrictions on the dfs
% are stricter. The inverse Wishart, for example, exists for nu > p-1, but
% its mean only exists for nu > p+1.
lb_df_barlett = (p-1);
lb_df_barlettL = (0:p-1);
lb_df_barlettU = flip(0:p-1);
lb_df_chi = 0;
x0_df_barlett = 2*p;
x0_df_barlettL = 2*p*ones(1,p);
x0_df_barlettU = 2*p*ones(1,p);
x0_df_chi = 5;

if strcmp( dist, 'Wish' )
    x0_df = x0_df_barlett;
    lb_df = lb_df_barlett;
elseif strcmp( dist, 'iWish' )
    x0_df = x0_df_barlett;  
    lb_df = lb_df_barlett;
elseif strcmp( dist, 'tWish' )
    x0_df = [ x0_df_barlett, x0_df_chi ];
    lb_df = [ lb_df_barlett, lb_df_chi ];
elseif strcmp( dist, 'itWish' )
    x0_df = [ x0_df_chi, x0_df_barlett ];
    lb_df = [ lb_df_chi, lb_df_barlett ];
elseif strcmp( dist, 'F' )
    x0_df = [ x0_df_barlett, x0_df_barlett ];   
    lb_df = [ lb_df_barlett, lb_df_barlett ];
elseif strcmp( dist, 'Riesz' )
    x0_df = x0_df_barlettL;
    lb_df = lb_df_barlettL;
elseif strcmp( dist, 'iRiesz2' )
    x0_df = x0_df_barlettU;
    lb_df = lb_df_barlettU;
elseif strcmp( dist, 'tRiesz' )
    x0_df = [ x0_df_barlettL, x0_df_chi ]; 
    lb_df = [ lb_df_barlettL, lb_df_chi ];
elseif strcmp( dist, 'itRiesz2' )
    x0_df = [ x0_df_chi, x0_df_barlettU ];
    lb_df = [ lb_df_chi, lb_df_barlettU];
elseif strcmp( dist, 'FRiesz' )
    x0_df = [ x0_df_barlettL, x0_df_barlettU ];
    lb_df = [ lb_df_barlettL, lb_df_barlettU];
elseif strcmp( dist, 'iFRiesz2' )
    x0_df = [ x0_df_barlettL, x0_df_barlettU ];
    lb_df = [ lb_df_barlettL, lb_df_barlettU ];
end
if isempty(x0)
    x0 = [vechchol(mean(R,3).*.25)', .05, .7, .01, .7, 2, x0_df];
end
                               % Stationarity
lb = [-inf(1,p_), zeros(1,4), 1, lb_df];
ub = [inf(1,p_), ones(1,4)];
Aeq = [zeros(1,p_+4), 1, zeros(1,length(x0_df))];
beq = 2;

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
            [],[],Aeq,beq, ...
            lb,ub,[], ...
            varargin{:} ...
        );
else
    [eparam,optimoutput] = ...
        my_fmincon(...
            @(param) obj_fun(param,1:p), ...
            x0', ...
            [],[],Aeq,beq, ...
            lb,ub,[], ...
            varargin{:} ...
        );
    optimoutput.perm_ = 1:p;
end
    
% Output Creation----------------------------------------------------------
[ nLogL, logLcontr, fit, fcst, eparam ] = obj_fun(eparam,optimoutput.perm_);
eparam.perm_ = optimoutput.perm_;

%[VCV,A,B,scores,hess,gross_scores] = robustvcv(fun, eparam, 3);
[VCV,scores,gross_scores] = vcv( obj_fun, eparam.all, optimoutput.perm_ );
tstats = eparam.all./sqrt(diag(VCV));


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

end