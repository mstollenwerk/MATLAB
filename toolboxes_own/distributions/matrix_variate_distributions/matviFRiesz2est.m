function [ eparam, tstats, logL, optimoutput ] = matviFRiesz2est( R, x0, varargin )
%
% Michael Stollenwerk
% michael.stollenwerk@live.com
% 28.04.2022
%
% DEPENDENCIES:
%
narginchk(2,inf);
[p,~,N] = size(R);
p_ = p*(p+1)/2;
%%
if ~isempty(varargin) && strcmp(varargin{1},'EstimationStrategy:MethodOfMoments')
    %% Optimization
    
    if isempty(x0)  
        x0 = [ 2*p.*ones(p,1); 2*p.*ones(p,1) ];
    end

    lb = [ flip((0:p-1)'); (0:p-1)' ];
    
    avgR = mean(R,3);

    if ~isempty(varargin(2)) && strcmp(varargin{2},'noperm')

        obj_fun = @(df) matvsiFRiesz2like( avgR, df(1:p), df(p+1:2*p), R );

        [eparam,optimoutput] = ...
            my_fmincon(...
                obj_fun, ...
                x0, ...
                [],[],[],[],lb,[],[],varargin{3:end} ...
            );

        obj_fun_opt_perm = obj_fun;  

        optimoutput.perm_ = 1:p;

    else

        obj_fun = @(df,perm_) matvsiFRiesz2like( avgR(perm_,perm_,:), ...
                                                 df(1:p), df(p+1:2*p), ...
                                                 R(perm_,perm_,:) );

        [eparam,optimoutput] = ...
            fmincon_Rieszperm(...
                p, ...
                obj_fun, ...
                x0, ...
                [],[],[],[],lb,[],[],varargin{2:end} ...
            );

        obj_fun_opt_perm = @(df) obj_fun(df,optimoutput.perm_);  
         
    end
    %% tstats
    %[VCV,A,B,scores,hess,gross_scores] = robustvcv(fun, eparam, 3);
    [VCV,scores,gross_scores] = vcv(obj_fun_opt_perm, eparam);

    tstats = eparam./sqrt(diag(VCV));

    tstats = struct(...
        'Sigma_', NaN(p_,1), ...
        'df_1', tstats(1 : p), ...               
        'df_2', tstats(p + 1 : p + p), ...               
        'all', [ NaN(p_,1); tstats ] ...
    );
    %% nLogL, logLcontr and eparam
    [nLogL, logLcontr, ~, ~, eparam ] = obj_fun_opt_perm( eparam );
    eparam.perm_ = optimoutput.perm_;

    aic = 2*nLogL + 2*(p_ + p + p);
    bic = 2*nLogL + log(N)*(p_ + p + p);
    logL = struct(...
        'logL', -nLogL, ...
        'logLcontr', logLcontr, ...
        'bic', bic, ...
        'aic', aic ...
    );
else
    %% Optimization
error("Not done yet!")     
%     error("Optimization over all parameters is done without optimization over asset ordering!")     
%     obj_fun = @(param) matviFRiesz2like( [], [], [], R, param );
% 
%     if isempty(x0)  
%         x0 = [ vech(chol(mean(R,3)/2/p, 'lower')); 2*p.*ones(p,1); 2*p.*ones(p,1) ];
%     end
% 
%     lb = [-inf(p_,1); (0:p-1)'; (p-(1:p)+2)']; %%%%%%%%%%%%
% 
%     [eparam,optimoutput] = ...
%         my_fmincon(...
%             obj_fun, ...
%             x0, ...
%             [],[],[],[],lb,[],[],varargin{:} ...
%         );
%     %% tstats
%     %[VCV,A,B,scores,hess,gross_scores] = robustvcv(fun, eparam, 3);
%     [VCV,scores,gross_scores] = ...
%         vcv(obj_fun, eparam);
% 
%     tstats = eparam./sqrt(diag(VCV));
% 
%     tstats = struct(...
%         'Sigma_', tstats(1:p_), ...
%         'df_1', tstats(p_ + 1 : p_ + p), ...               
%         'df_2', tstats(p_ + p + 1 : p_ + p + p), ...               
%         'all', tstats ...
%     );
%     %% nLogL, logLcontr and eparam
%     [nLogL, logLcontr, ~, ~, eparam ] = obj_fun( eparam );
% 
%     aic = 2*nLogL + 2*(p_ + p + p);
%     bic = 2*nLogL + log(N)*(p_ + p + p);
%     logL = struct(...
%         'logL', -nLogL, ...
%         'logLcontr', logLcontr, ...
%         'bic', bic, ...
%         'aic', aic ...
%     );
end

end
