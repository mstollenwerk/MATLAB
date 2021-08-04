function [eparam, tstats, logL, optimoutput] = matvncFBOPSest(X, mhg_precision, x0, varargin)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

%% Input Checking
% Will be added later
[p,~,N] = size(X);
p_ = p*(p+1)/2;
%% Optimization
if isempty(x0)
    %"Pre-Optimization" based on central F.
    x0 = matvFest(X, [], 'Display', 'off', varargin{:});
    %"Pre-Optimization" to get good "Sigma_ to Omega_ ratio".
    meanX = x0.Sigma_*x0.df_1/(x0.df_2 - 2);
    weight = .25:.05:.85;
    y = NaN(length(weight),1);
    for ii = 1:length(weight) 
        y(ii) = matvncFBOPSlike( ...
            meanX*weight(ii), ...
            meanX*(1 - weight(ii))*x0.df_1, ...
            x0.df_1, x0.df_2 + p - 1, mhg_precision, X ...
        ); 
    end
    
    y = y.*~isinf(y);
    [~, min_ii] = min(y);
    x0 = [vechchol(meanX*weight(min_ii));
          vechchol(meanX*(1-weight(min_ii))*x0.df_1);
          x0.df_1;
          x0.df_2 + p - 1];
      
% You should start far away from Omega_ = 0-Matrix!
%     x0 = [vechchol(eye(p).*.01);vechchol(mean(X,3)*2*p);2*p;2*p];
%     x0 = [vechchol(mean(X,3));vechchol(eye(p).*.1);2*p;2*p];
end

logL_improvement = inf;
ii = 1;
[eparam,optimoutput{1}] = ...
    my_fmincon(...
        @(param) matvncFBOPSlike([],[],[],[],mhg_precision,X,param),...
        x0,...        
        [],[],[],[],[-inf(2*p_,1);p-1;p-1],[],[],...
        varargin{:}...
    );
nlogL = optimoutput{1}.history.fval(end);
while logL_improvement > 1e-6
    
    disp('--------------------------------------')
    disp('--------------------------------------')
    disp('--------------------------------------')
    disp(['Beginning partial optimizations # ' num2str(ii) '.'])
    disp('--------------------------------------')
    disp('--------------------------------------')
    disp('--------------------------------------')
    
    free_indices = ...
        { [1 : p_, p_ + p_ + 1, p_ + p_ + 2], ... % Sigma, DOFs
          [p_ + 1 : p_ + p_, p_ + p_ + 1, p_ + p_ + 2] }; % Omega, DOFs

    [eparam_new,optimoutput{2,ii}] = ...
        my_fminunc_partial( ...
            @(param) matvncFBOPSlike([],[],[],[],mhg_precision,X,param), ...
            eparam,free_indices,varargin{:} ...
        );
    
    disp('--------------------------------------')
    disp('--------------------------------------')
    disp('--------------------------------------')
    disp(['Beginning full optimization # ' num2str(ii+1) '.'])
    disp('--------------------------------------')
    disp('--------------------------------------')
    disp('--------------------------------------')
    [eparam_new,optimoutput{1,ii+1}] = ...
        my_fmincon(...
                @(param) matvncFBOPSlike([],[],[],[],mhg_precision,X,param),...
                eparam_new,...        
                [],[],[],[],[-inf(2*p_,1);p-1;p-1],[],[],...
                varargin{:}...
            ); 
    
    logL_improvement = nlogL - optimoutput{1,ii+1}.history.fval(end);
    if logL_improvement > 1e-6
        nlogL = optimoutput{1,ii+1}.history.fval(end);
        eparam = eparam_new;
    else
        optimoutput = optimoutput(1:end-3);
    end
    ii = ii + 1;
end
%% tstats
%[VCV,A,B,scores,hess,gross_scores] = robustvcv(fun, eparam, 3);
[VCV,scores,gross_scores] = ...
    vcv(@(param) matvncFBOPSlike([],[],[],[],mhg_precision,X,param), eparam);

tstats = eparam./sqrt(diag(VCV));

tstats = struct(...
    'Sigma_', tstats(1:p_), ...
    'Omega_', tstats(p_+1:p_+p_), ...
    'df_1', tstats(p_ + p_ + 1), ...           
    'df_2', tstats(p_ + p_ + 2), ...     
    'all', tstats ...
);
%% nLogL, logLcontr and eparam  
[nLogL, logLcontr, ~, ~, eparam] = matvncFBOPSlike([],[],[],[],mhg_precision,X,eparam);

aic = 2*nLogL + 2*numel(x0);
bic = 2*nLogL + log(N)*numel(x0);
logL = struct(...
    'logL', -nLogL, ...
    'logLcontr', logLcontr, ...
    'bic', bic, ...
    'aic', aic ...
);
end

