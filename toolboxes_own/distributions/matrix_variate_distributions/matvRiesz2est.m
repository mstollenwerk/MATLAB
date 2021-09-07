function [ eparam, tstats, logL, optimoutput ] = matvRiesz2est( data_mat, x0, varargin )
%MATVRIESZEST
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
% 11.02.2021
%
% DEPENDENCIES:
%
narginchk(2,inf);
[p,~,N] = size(data_mat);
p_ = p*(p+1)/2;
%%
if ~isempty(varargin) && strcmp(varargin{1},'EstimationStrategy:MethodOfMoments')
    %% Optimization
    
    obj_fun = @(df) matvRiesz2like( cholU(mean(data_mat,3))/diag(df)*cholU(mean(data_mat,3))', df, data_mat );

    if isempty(x0)  
        x0 = 2*p.*ones(p,1);
    end

    lb = [fliplr(0:p-1)]';

    [eparam,optimoutput] = ...
        my_fmincon(...
            obj_fun, ...
            x0, ...
            [],[],[],[],lb,[],[],varargin{2:end} ...
        );
    optimoutput.perm_ = 1:p;
    
    perm_improvement = inf;
    useless_counter = 0;
    while perm_improvement > 0
        useless_counter = useless_counter + 1;
        rng(useless_counter) % For comparison, always try the same permutations, but not the same in every while loop iteration.
        for ii = 1:100
            perm_(ii,:) = randperm(p);
        end
        perm_ = unique(perm_, 'rows', 'stable');
        perm_ = setdiff(perm_, optimoutput.perm_, 'rows'); %remove previously optimized 
        for ii = 1:size(perm_,1)
            nLogL_permuted_assets(ii) = ...
                matvRiesz2like( ...
                    cholU(mean(data_mat(perm_(ii,:),perm_(ii,:),:),3))/diag(eparam)*cholU(mean(data_mat(perm_(ii,:),perm_(ii,:),:),3))', ...
                    eparam, data_mat(perm_(ii,:),perm_(ii,:),:) ...
                );
        end
        [~,min_ii] = min(nLogL_permuted_assets);
        
        disp(strcat("Optimizing over asset permutation (",num2str(perm_(min_ii,:)),")"))
        obj_fun_min_ii = @(df) matvRiesz2like( ...
            cholU(mean(data_mat(perm_(min_ii,:),perm_(min_ii,:),:),3))/diag(df)*cholU(mean(data_mat(perm_(min_ii,:),perm_(min_ii,:),:),3))', ...
            df, data_mat(perm_(min_ii,:),perm_(min_ii,:),:) );
        [eparam_min_ii,optimoutput_min_ii] = ...
            my_fmincon(...
                obj_fun_min_ii, ...
                eparam, ...
                [],[],[],[],lb,[],[],varargin{2:end} ...
            );
        
        perm_improvement = -optimoutput_min_ii.history.fval(end) + optimoutput.history.fval(end);
        if perm_improvement > 0
            perm_improvement
            obj_fun = obj_fun_min_ii;
            eparam = eparam_min_ii;
            optimoutput = optimoutput_min_ii;
            optimoutput.perm_ = perm_(min_ii,:);
        else
            disp("No likelihood improvement with new asset permuatation.")
        end
        clear perm_ nLogL_permuted_assets
    end  
    %% tstats
    %[VCV,A,B,scores,hess,gross_scores] = robustvcv(fun, eparam, 3);
    [VCV,scores,gross_scores] = vcv(obj_fun, eparam);

    tstats = eparam./sqrt(diag(VCV));

    tstats = struct(...
        'Sigma_', NaN(p_,1), ...
        'df', tstats, ...               
        'all', [NaN(p_,1); tstats] ...
    );
    %% nLogL, logLcontr and eparam
    [nLogL, logLcontr, ~, ~, eparam ] = obj_fun( eparam );

    aic = 2*nLogL + 2*(p_ + p);
    bic = 2*nLogL + log(N)*(p_ + p);
    logL = struct(...
        'logL', -nLogL, ...
        'logLcontr', logLcontr, ...
        'bic', bic, ...
        'aic', aic ...
    );
else
    %% Optimization
    warning("Optimization over all parameters is done without optimization over asset ordering!")
    obj_fun = @(param) matvRiesz2like( [], [], data_mat, param );

    if isempty(x0)  
        x0 = [ vech(chol(mean(data_mat,3)/2/p, 'lower')); 2*p.*ones(p,1) ];
    end

    lb = [-inf(p_,1);flipud((0:p-1)')];

    [eparam,optimoutput] = ...
        my_fmincon(...
            obj_fun, ...
            x0, ...
            [],[],[],[],lb,[],[],varargin{:} ...
        );
    %% tstats
    %[VCV,A,B,scores,hess,gross_scores] = robustvcv(fun, eparam, 3);
    [VCV,scores,gross_scores] = vcv(obj_fun, eparam);

    tstats = eparam./sqrt(diag(VCV));

    tstats = struct(...
        'Sigma_', tstats(1:p_), ...
        'df', tstats(p_ + 1 : p_ + p), ...               
        'all', tstats ...
    );
    %% nLogL, logLcontr and eparam
    [nLogL, logLcontr, ~, ~, eparam ] = obj_fun( eparam );

    aic = 2*nLogL + 2*(p_ + p);
    bic = 2*nLogL + log(N)*(p_ + p);
    logL = struct(...
        'logL', -nLogL, ...
        'logLcontr', logLcontr, ...
        'bic', bic, ...
        'aic', aic ...
    );
end
end
