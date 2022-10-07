function [eparam, tstats, logL, fit, fcst, optimoutput] = ...
	caicov_scalar_BEKK_estim_targeting_new( R, p, q, dist, x0, varargin )
% THIS DOES NOT WORK BECAUSE:
% You need the gradient and hessian of the STANDARDIZED
% distributoins. These I cannot get for the FRiesz.
%
% Michael Stollenwerk
% michael.stollenwerk@live.com
% 05.04.2021

%% Input Checking
% Will be added later
[k,~,T] = size(R);
k_ = k*(k+1)/2;

%% Optimization
obj_fun = @(param, perm_) obj_fun_wrapper(param, R(perm_,perm_,:), p, q, dist);
% x0-----------------------------------------------------------------------
if ~isempty(x0) && (isinf(obj_fun(x0,1:k)) || isnan(obj_fun(x0,1:k)))
    error('User supplied x0 invalid.')
end
% The lower bounds below are for existence of the distributions. For
% existence of the expected value matrix Sigma_, restrictions on the dfs
% are stricter. The inverse Wishart, for example, exists for nu > p-1, but
% its mean only exists for nu > p+1.
lb_df_barlett = [(k-1)];
lb_df_barlettL = [(0:k-1)];
lb_df_barlettU = [flip(0:k-1)];
lb_df_chi = [0];
x0_df_barlett = [k+10];
x0_df_barlettL = [k*ones(1,k)+10];
x0_df_barlettU = [k*ones(1,k)+10];
x0_df_chi = [5];

if strcmp( dist, 'Wish' )
    x0_df = x0_df_barlett;
    lb_df = lb_df_barlett;
%     Aeq = [ zeros(4,6) eye(4) ];
%     beq = zeros(4,1);       
elseif strcmp( dist, 'iWish' )
    x0_df = x0_df_barlett;  
    lb_df = lb_df_barlett;
%     Aeq = [ zeros(4,5*p+1) eye(4) ];
%     beq = zeros(4,1);    
elseif strcmp( dist, 'tWish' )
    x0_df = [ x0_df_barlett, x0_df_chi ];
    lb_df = [ lb_df_barlett, lb_df_chi ];
%     Aeq = [ zeros(4,6) eye(4) zeros(4,5) ;
%             zeros(4,11) eye(4)];
%     beq = zeros(8,1);    
elseif strcmp( dist, 'itWish' )
    x0_df = [ x0_df_chi, x0_df_barlett ];
    lb_df = [ lb_df_chi, lb_df_barlett ];
%     Aeq = [ zeros(4,6) eye(4) zeros(4,5) ;
%             zeros(4,11) eye(4)];
%     beq = zeros(8,1);      
elseif strcmp( dist, 'F' )
    x0_df = [ x0_df_barlett, x0_df_barlett ];   
    lb_df = [ lb_df_barlett, lb_df_barlett ]; 
%     Aeq = [ zeros(4,6) eye(4) zeros(4,5) ;
%             zeros(4,11) eye(4)];
%     beq = zeros(8,1);        
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
    x0 = [0.05, 0.9, x0_df]';
end
% A = [ ones(1,p+q) zeros(1, length(x0)-p-q) ];   % Stationarity
% b = 1;                                          % Stationarity
lb = [0, 0, lb_df];

if contains(dist,'Riesz')
    perm_optim = 1;
else
    perm_optim = 0;
end
if ~isempty(varargin) && strcmp(varargin{1},'noperm')
    perm_optim = 0;
    varargin = varargin(2:end);
end
varargin = [varargin, {'Algorithm','trust-region','SpecifyObjectiveGradient', true, 'HessianFcn', 'objective'}];
if perm_optim
    [eparam,optimoutput] = ...
        fminunc_Rieszperm(...
            k, ...
            obj_fun, ...
            x0, ...
            varargin{:} ...
        );
else
    [eparam,optimoutput] = ...
        my_fminunc(...
            @(param) obj_fun(param,1:k), ...
            x0, ...
            varargin{:} ...
        );
    optimoutput.perm_ = 1:k;
end

% Output Creation----------------------------------------------------------
[ nLogL, logLcontr, g, ~, Sigma_, eparam ] = obj_fun( eparam, optimoutput.perm_ );
eparam.perm_ = optimoutput.perm_;

% Tstats-------------------------------------------------------------------
B=cov(g);%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
VCV=inv(B)/T;
tstats = eparam.all(k_+1:end)./sqrt(diag(VCV));

aic = 2*nLogL + 2*(numel(x0)+k_);
bic = 2*nLogL + log(T)*(numel(x0)+k_); % see Yu, Li and Ng (2017)
logL_detR_part = -(k+1)/2*sum(logdet3d(R));
logL = struct(...
    'logL', -nLogL, ...
    'logL_detR_part', logL_detR_part, ...    
    'aic', aic, ...
    'bic', bic, ...
    'logLcontr', logLcontr ...
);

fit = struct( ...
    'Sigma_', Sigma_(:,:,1:T) ...
);
fcst = struct('Sigma_', Sigma_(:,:,T+1:end));

end

function [ nLogL, g, H, logLcontr, Sigma_, param_out ] = obj_fun_wrapper(param, R, p, q, dist) 
    
    T = size(R,3);

    param_g_arch = param(1:p+q);
    param_df = param(p+q+1:end);
    
    if sum(param_g_arch) >= 1
        nLogL = inf;
        g = NaN;
        H = NaN;
        logLcontr = NaN;
        Sigma_ = NaN;
        param_out = NaN;
        return
    end
    
    meanR = mean(R,3);
    
    [ nLogL, logLcontr, Sigma_, param_out ] = caicov_scalar_BEKK_likeRec( ...
        [vechchol( meanR*(1-sum(param(1:p+q))) ); param], ...
        p, ...
        q, ...
        R, ...
        dist ...
    );

    if isinf(nLogL)
        g = NaN;
        H = NaN;
        return
    end

    Omega_ = matvStandardize(dist, Sigma_(:,:,1:T), param_df);
%     g_df = matvLogLikeGradientTheta('sRiesz', Sigma_(:,:,1:T), repmat(param_df',T,1), R);
    g_df = matvLogLikeGradientTheta(dist, Omega_, repmat(param_df',T,1), R);
        
    H_df = matvLogLikeHessianTheta(dist, Omega_, repmat(param_df',T,1), R);
    
    [g_num, H_g_garch_ALL] = grad_hessian_2sided_ij(...
        @(x) caicov_scalar_BEKK_likeRec( ...
                [vechchol( meanR*(1-sum(x(1:p+q))) ); x],p,q,R,dist ...
             ), ...
        param, ...
        [1:p+q], ...
        [1:length(param)] ...
    );

    g = [g_num(1:p+q); -sum(g_df,1)'];

    H = [ H_g_garch_ALL(1:p+q,:);
          H_g_garch_ALL(1:p+q,p+q+1:end)', -sum(H_df,3)];
%           H_df_g_garch(p+q+1:end,1:p+q), sum(H_df,3)];
    
end