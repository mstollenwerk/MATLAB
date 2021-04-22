function [eparam, tstats, logL, fit, fcst, optimoutput] = ...
	gas_ncF_estim_staticOmega( data_RC, p, q, x0, mhg_precision, varargin )
%% Input Checking
% Will be added later
[k,~,TT] = size(data_RC);
if k > 5
    warning('For k > 5 estimation might be too slow due to curse of dimensionality.')
end
%% Optimization
% Objective Function-------------------------------------------------------
fun = @( param ) gas_ncF_likeRec_staticOmega( param, p, q, data_RC, mhg_precision );
% x0-----------------------------------------------------------------------
if isempty(x0) || isstruct(x0)
    if isfield(x0,'gas_intrcpt')
        gas_intrcpt = x0.gas_intrcpt;
    else
        gas_intrcpt = mean(data_RC,3).*.1;
    end
    if isfield(x0,'gas_arch_matrix')
        gas_arch_matrix = x0.gas_arch_matrix;
    else
        gas_arch_matrix = repmat(zeros(k),1,1,p);
    end
    if isfield(x0,'gas_garch_matrix')
        gas_garch_matrix = x0.gas_garch_matrix;
    else
        gas_garch_matrix = repmat(zeros(k),1,1,q);
    end
    if isfield(x0,'Omega_')
        Omega_ = x0.Omega_;
    else
        Omega_ = mean(data_RC,3)*2*k;        
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
    x0 = [ vechchol(gas_intrcpt);
           gas_arch_matrix(:);               
           gas_garch_matrix(:);              
           vechchol(Omega_);                    
           df_F_1;                           
           df_F_2;                                                                
    ];
else
    error('x0 must be empty or struct.')
end
% Optimization-------------------------------------------------------------
[eparam,optimoutput,] = my_fminunc( fun,x0,varargin{:} );
% Tstats-------------------------------------------------------------------
%[VCV,A,B,scores,hess,gross_scores] = robustvcv(fun, eparam, 3);
[VCV,scores,gross_scores] = vcv(fun, eparam);
tstats = eparam./sqrt(diag(VCV));

% Output Creation----------------------------------------------------------

[ nLogL, logLcontr, Sigma_, ScoreMat, eparam ] = fun(eparam);

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
    'eExp', Sigma_(:,:,1:end-1) + eparam.Omega_/eparam.df_F_1,...
    'ScoreMat', ScoreMat(:,:,1:end)...
);
fcst = struct(...
    'Sigma_', Sigma_(:,:,end),...
    'eExp', Sigma_(:,:,end) + eparam.Omega_/eparam.df_F_1 ...
);
end