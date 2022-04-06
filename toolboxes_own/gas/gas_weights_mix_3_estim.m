function [eparam, tstats, logL, fit, fcst, weights, optimoutput] = ...
	gas_weights_mix_3_estim( logL1, logL2, logL3, x0, varargin )
%
%
% Michael Stollenwerk
% michael.stollenwerk@live.com
% 05.04.2022

%% Optimization
obj_fun = @(param) gas_scalar_BEKK_mix_titF_Riesz_likeRec( param, logL1, logL2, logL3 );
% x0-----------------------------------------------------------------------
if ~isempty(x0) && (isinf(obj_fun(x0,1:k)) || isnan(obj_fun(x0,1:k)))
    error('User supplied x0 invalid.')
end
if isempty(x0)
    x0 = repmat([.33 .33 .95 .0001]',2,1);
end
% Restrictions-------------------------------------------------------------
lb = [0 0 0 -inf 0 0 0 -inf]';
ub = [1 1 inf inf 1 1 inf inf]';
% Optimization-------------------------------------------------------------
options_ = optimoptions('fmincon',varargin{:});
[eparam,optimoutput] = ...
    fmincon(...
        obj_fun, ...
        x0, ...
        [],[],[],[], ...
        lb,ub,[], ...
        options_ ...
    );

% Output Creation----------------------------------------------------------
[ nLogL, logLcontr, weights, scores_weights ] = obj_fun( eparam );

%[VCV,A,B,scores,hess,gross_scores] = robustvcv(fun, eparam, 3);
[VCV,scores,gross_scores] = vcv( obj_fun, eparam );
tstats = eparam./sqrt(diag(VCV));


% aic = 2*nLogL + 2*(numel(x0)+k_);
% bic = 2*nLogL + log(T)*(numel(x0)+k_); % see Yu, Li and Ng (2017) 
% logL_detR_part = -(k+1)/2*sum(logdet3d(R));
logL = struct(...
    'logL', -nLogL,... 'logL_detR_part', logL_detR_part, ...'aic', aic,...'bic', bic,...
    'logLcontr', logLcontr...
);

T = length(logL1);
weights_fit{1} = weights{1}(1:T);
weights_fit{2} = weights{2}(1:T);
fit = struct( ...
    'weights', weights_fit, ...
    'score_weights', scores_weights ...
);

weights_fcst{1} = weights{1}(T+1:end);
weights_fcst{2} = weights{2}(T+1:end);
fcst = struct('weights', weights_fcst);

end