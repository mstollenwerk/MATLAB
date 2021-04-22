function [eparam, tstats, logL, fit, fcst, optimoutput] = ...
	gas_ncF_estim( data_RC, p, q, r, x0, mhg_precision, varargin )
%% Input Checking
% Will be added later
[k,~,TT] = size(data_RC);
k_ = k*(k+1)/2;
if k > 5
    warning('For k > 5 estimation might be too slow due to curse of dimensionality.')
end
%% Optimization
% Objective Function-------------------------------------------------------
fun = @( param ) gas_ncF_likeRec( param, p, q, r, data_RC, mhg_precision );
% x0-----------------------------------------------------------------------
if isempty(x0) || isstruct(x0)
    if isfield(x0,'gas_intrcpt')
        gas_intrcpt = x0.gas_intrcpt;
    else
        gas_intrcpt = mean(data_RC,3)*.1;
    end
    if isfield(x0,'gas_arch_matrix')
        gas_arch_matrix = x0.gas_arch_matrix;
    else
        gas_arch_matrix = repmat(eye(k)*.1/p,1,1,p);
    end
    if isfield(x0,'gas_garch_matrix')
        gas_garch_matrix = x0.gas_garch_matrix;
    else
        gas_garch_matrix = repmat(eye(k)*.5/q,1,1,q);
    end
    if isfield(x0,'war_matrix')
        war_matrix = x0.war_matrix;
    else
        war_matrix = repmat(eye(k)*.01,1,1,r);        
    end
    if isfield(x0,'df_F_1')
        df_F_1 = x0.df_F_1;
    else
        df_F_1 = 2*k;
    end    
    if isfield(x0,'df_F_2')
        df_F_2 = x0.df_F_2;
    else
        df_F_2 = 2*k;
    end 
    x0 = [ vech(chol(gas_intrcpt, 'lower'));
           gas_arch_matrix(:);               
           gas_garch_matrix(:);              
           war_matrix(:);                    
           df_F_1;                           
           df_F_2;                                                                
    ];
else
    error('x0 must be empty or struct.')
end
% Restrictions-------------------------------------------------------------
lb = [ -inf(k_,1);
       -inf((p+q+r)*k^2,1)
% Restrictions below are for uniqueness of (G)ARCH Matrices. (not
% necessary)
%        kron(ones(p,1),[0;-inf(k^2-1,1)]);
%        kron(ones(q,1),[0;-inf(k^2-1,1)]);	   
%        kron(ones(r,1),[0;-inf(k^2-1,1)]);
       k+eps;k+1+eps;2+eps
];
% Optimization-------------------------------------------------------------
[eparam,optimoutput,] = my_fmincon(...
    fun,x0,[],[],[],[],lb,[],[],varargin{:});
% Tstats-------------------------------------------------------------------
%[VCV,A,B,scores,hess,gross_scores] = robustvcv(fun, eparam, 3);
[VCV,scores,gross_scores] = vcv(fun, eparam);
tstats = eparam./sqrt(diag(VCV));

% Output Creation----------------------------------------------------------

[ nLogL, logLcontr, Sigma_, Omega_, ScoreMat, eparam ] = fun(eparam);

aic = 2*nLogL + 2*(numel(x0));
bic = 2*nLogL + log(TT)*numel(x0); % see Yu, Li and Ng (2017) 
logL = struct(...
    'logL', -nLogL,...
    'aic', aic,...
    'bic', bic,...
    'logLcontr', logLcontr...
);
fit = struct(...
    'Sigma_', Sigma_(:,:,1:end-1),...
    'Omega_', Omega_(:,:,1:end-1),...
    'eExp', Sigma_(:,:,1:end-1) + Omega_(:,:,1:end-1)/eparam.df_F_1,...
    'ScoreMat', ScoreMat(:,:,1:end)...
);
fcst = struct(...
    'Sigma_', Sigma_(:,:,end),...
    'Omega_', Omega_(:,:,end),...
    'eExp', Sigma_(:,:,end) + Omega_(:,:,end)/eparam.df_F_1 ...
);
end