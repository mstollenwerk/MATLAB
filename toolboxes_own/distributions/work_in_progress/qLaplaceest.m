function [ eparam, tstats, logL, optimoutput ] = matvLaplacWishest( X, x0, varargin )
%MATVLAPLACEWISHEST
%
% USAGE:
%   
%
% INPUTS:
%   
%
% OUTPUTS:
%   
%  See also 
%
% COMMENTS:
%   
% REFERENCES:
%      [1] Stollenwerk (2020)
%
% Michael Stollenwerk
% michael.stollenwerk@live.com
% 18.02.2020
%
% DEPENDENCIES:
%
%% Input Checking
% Will be added later
[p,~,N] = size(X);
p_ = p*(p+1)/2;
if ~isempty(varargin) && strcmp(varargin{1},'EstimationStrategy:MethodOfMoments')
    %% Optimization
    obj_fun = @(df) matvLaplaceWishlike( mean(X,3)*(df(1)-p-1)/df(2),df(1),df(2),X);

    if isempty(x0)
        x0 = [2*p; 2*p];
    end

    [eparam,optimoutput] = ...
        my_fmincon(...
            obj_fun,...
            x0,...
            [],[],[],[],[p+1;p-1],[],[],...
            varargin{2:end}...
        );
    %% tstats
    %[VCV,A,B,scores,hess,gross_scores] = robustvcv(fun, eparam, 3);
    [VCV,scores,gross_scores] = vcv(obj_fun, eparam);

    tstats = eparam./sqrt(diag(VCV));

    tstats = struct(...
        'Sigma_', NaN(p_,1), ...
        'df_1', tstats(1), ...           
        'df_2', tstats(2), ...     
        'all', [NaN(p_,1); tstats] ...
    );
    %% nLogL, logLcontr and eparam  
    [nLogL, logLcontr, ~, ~, eparam] = obj_fun(eparam);

    aic = 2*nLogL + 2*(p_ + 2);
    bic = 2*nLogL + log(N)*(p_ + 2);
    logL = struct(...
        'logL', -nLogL, ...
        'logLcontr', logLcontr, ...
        'bic', bic, ...
        'aic', aic ...
    );
else
    %% Optimization
    obj_fun = @(param) matvLaplaceWishlike([],[],[],X,param);

    if isempty(x0)
        x0 = [vechchol( mean(X,3)*(p-1)/p ); 2*p; 2*p];
    end

    [eparam,optimoutput] = ...
        my_fmincon(...
            obj_fun,...
            x0,...
            [],[],[],[],[-inf(p_,1);p+1;p-1],[],[],...
            varargin{:}...
        );
    %% tstats
    %[VCV,A,B,scores,hess,gross_scores] = robustvcv(fun, eparam, 3);
    [VCV,scores,gross_scores] = vcv(obj_fun, eparam);

    tstats = eparam./sqrt(diag(VCV));

    tstats = struct(...
        'Sigma_', tstats(1:p_), ...
        'df_1', tstats(p_ + 1), ...           
        'df_2', tstats(p_ + 2), ...     
        'all', tstats ...
    );
    %% nLogL, logLcontr and eparam  
    [nLogL, logLcontr, ~, ~, eparam] = obj_fun(eparam);

    aic = 2*nLogL + 2*(p_ + 2);
    bic = 2*nLogL + log(N)*(p_ + 2);
    logL = struct(...
        'logL', -nLogL, ...
        'logLcontr', logLcontr, ...
        'bic', bic, ...
        'aic', aic ...
    );
end
end
% 
% %% Error checking
% narginchk(1,inf);
% if nargin==1
%     x0 = [];
% else
%     x0 = varargin{1};
% end
% [p,p2,N] = size(X);
% if p~=p2
%     error('Covariance matrix data must be square.')
% end
% p_ = p*(p+1)/2;
% %% Optimization
% % Objective Function ------------------------------------------------------
% fun = @( param ) like_wrapper(param, X, p_, N);
% % x0 ----------------------------------------------------------------------
% if isempty(x0) || isstruct(x0) 
%     if isfield(x0,'df')
%         df = x0.df;
%     else
%         df = 2*p;
%     end    
%     if isfield(x0,'Sigma_')
%         Sigma_ = x0.Sigma_;
%     else
%         Sigma_ = mean(X,3)*df;
%     end      
%     x0 = [ vech(chol(Sigma_, 'lower')); df ]'; % k*(k+1)/2 + 1
%     try
%         fun(x0);
%     catch
%         error('x0 is not valid.')
%     end
% else
%     error('x0 must be empty or struct.')
% end
% % Restrictions-------------------------------------------------------------
% % Optimoptions-------------------------------------------------------------
% options = optimoptions('fmincon',...
%     'MaxFunctionEvaluations', 1e5 ...
% );
% options = optimoptions(options, varargin{2:end});
% % Restrictions-------------------------------------------------------------
% A = []; %
% b = [];
% Aeq = [];
% beq = [];
% lb = [-inf(p_,1);p-1];
% ub = [];
% nonlcon = [];
% % Optimization-------------------------------------------------------------
% warning('off','MATLAB:nearlySingularMatrix') %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% tic
% [eparam,~,exitflag,optimoutput,grad,hessian] = fmincon(...
%     fun,x0,A,b,Aeq,beq,lb,ub,nonlcon,options);
% optimoutput.estimation_time = toc;
% optimoutput.exitflag = exitflag;
% optimoutput.grad = grad;
% optimoutput.hessian = hessian;
% warning('on','MATLAB:nearlySingularMatrix') %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Tstats-------------------------------------------------------------------
% %[VCV,A,B,scores,hess,gross_scores] = robustvcv(fun, eparam, 3);
% [VCV,scores,gross_scores] = vcv(fun, eparam);
% tstats = eparam./sqrt(diag(VCV)');
% % Output Creation----------------------------------------------------------
% [nLogL, logLcontr, Sigma_, df ] = like_wrapper(eparam, X, p_, N);
% 
% eparam = struct( ...
%     'Sigma_', Sigma_, ...
%     'df', df, ...          
%     'all', eparam ...
% );
% 
% aic = 2*nLogL + 2*(numel(x0));
% bic = 2*nLogL + log(N)*(numel(x0));
% logL = struct(...
%     'logL', -nLogL, ...
%     'logLcontr', logLcontr, ...
%     'bic', bic, ...
%     'aic', aic ...
% );
% tstats = struct(...
%     'Sigma_', tstats(1:p*(p+1)/2),...
%     'df', tstats(p*(p+1)/2 + 1) ...           
% );
% %% Likelihood wrapper function
% function [nLogL, logLcontr, Sigma_, df ] = like_wrapper(param, data_mat, p_, N) 
% 
% Sigma_ = ivech(param(1:p_), 'lower');
% Sigma_ = Sigma_*Sigma_';
% df = param(p_+1);
% 
% [ nLogL, logLcontr ] = qLaplacelike( ...
%     data_mat, ...
%     ones(N,1)*df, ...
%     repmat(Sigma_,1,1,N) ...
% );
% 
% end
% end