function [eparam, tstats, logL, fit_fcst, optimoutput] = ...
	gas_estim( dist, data_RC, x0, varargin )
%
% Michael Stollenwerk
% michael.stollenwerk@live.com
% 11.08.2020

%% Input Checking
% Will be added later
[p,~,N] = size(data_RC);
p_ = p*(p+1)/2;
%% Optimization
if strcmp( dist, 'Wish' )
    if isempty(x0) 
        x0 = [ vech(chol( mean(data_RC,3)*(1-0.1) ,'lower'));
               0.1;
               0.7;
               2*p
        ];
    end
    A = [ zeros(1, p_), 1, -1, 0 ];
    lb = [-inf(1, p_), -eps, -eps, p-1];  
    ub = [ inf(1, p_), 1, 1, inf ];
elseif strcmp( dist, 'iWish' )
    if isempty(x0) 
        x0 = [ vech(chol( mean(data_RC,3)*(1-0.1) ,'lower'));
               0.01;
               0.7;
               2*p
        ];
    end
    A = [ zeros(1, p_), 1, -1, 0 ];
    lb = [-inf(1, p_), -eps, -eps, p+1];  
    ub = [ inf(1, p_), 1, 1, inf ];     
elseif strcmp( dist, 'F' )
    if isempty(x0) 
        x0 = [ vech(chol( mean(data_RC,3)*(1-0.1) ,'lower'));
               0.1;
               0.7;
               2*p;
               p
        ];
    end
    A = [ zeros(1, p_), 1, -1, 0, 0 ];
    lb = [-inf(1, p_), -eps, -eps, p-1, 0];
    ub = [ inf(1, p_), 1, 1, inf, inf ];
elseif strcmp( dist, 'qt' )
    if isempty(x0) 
        x0 = [ vech(chol( mean(data_RC,3)*(1-0.1) ,'lower'));
               0.1;
               0.7;
               10;
               2*p
        ];
    end
    A = [ zeros(1, p_), 1, -1, 0, 0 ];
    lb = [-inf(1, p_), -eps, -eps, p-1, 0];
    ub = [ inf(1, p_), 1, 1, inf, inf ];
elseif strcmp( dist, 'qtGammaprod' )
    if isempty(x0) 
        x0 = [ vech(chol( mean(data_RC,3)*(1-0.1) ,'lower'));
               0.1;
               0.7;
               10;
               2*p;
               10
        ];
    end
    A = [ zeros(1, p_), 1, -1, 0, 0, 0 ];
    lb = [-inf(1, p_), -eps, -eps, p-1, 0, 0];
    ub = [ inf(1, p_), 1, 1, inf, inf, inf ];          
end


[eparam,optimoutput] = ...
	my_fmincon(...
		@(param) gas_likeRec( param, data_RC, dist), ... % objective function
		x0, ...                            % x0
        A, ...                             % A: score_param <= garch_param
        -eps, ...                          % b
        [],[], ...                         % Aeq, Beq 
        lb, ...                            % lb
        ub, ...                            % ub
        [], ...                            % nonlcon
		varargin{:} ...                    % options
	);
%% tstats
%[VCV,A,B,scores,hess,gross_scores] = robustvcv(fun, eparam, 3);
[VCV,scores,gross_scores] = ...
	vcv( ...
		@(param) gas_likeRec( param, data_RC, dist), ...
		eparam ...
	);

tstats = eparam./sqrt(diag(VCV));
%% eparam, logL, fit_fcst 
[ nLogL, logLcontr, Sigma_, Score, eparam ] = ...
    gas_likeRec( eparam, data_RC, dist);

aic = 2*nLogL + 2*numel(x0);
bic = 2*nLogL + log(N)*numel(x0); % see Yu, Li and Ng (2017) 
logL = struct(...
    'logL', -nLogL,...
    'aic', aic,...
    'bic', bic,...
    'logLcontr', logLcontr...
);

fit_fcst = struct('Sigma_', Sigma_, 'Score', Score);
end