function [eparam, logL, fit, fcst, optimoutput] = cancF_estim(data_RC, x0)

%% Input Checking
% Will be added later
[k,~,T] = size(data_RC);
%% Optimization
fun = @( param ) likeRec_wrapper(param, data_RC);
if isempty(x0)
    x0 = [ vech(eye(k)); % intrcpt       k*(k+1)/2
           .1;           % arch_param    k*(k+1)/2 + 1
           .8;           % garch_param   k*(k+1)/2 + 1 + 1
           vech(eye(k))  % war_matrix    k*(k+1)/2 + 1 + 1 + k*(k+1)/2
           2*k;          % df_F_1        k*(k+1)/2 + 1 + 1 + k*(k+1)/2 + 1
           2*k           % df_F_2        k*(k+1)/2 + 1 + 1 + k*(k+1)/2 + 2
    ];
elseif isstruct(x0)
    disp('Converting struct x0 to vector')
    x0 = [ vech(chol(x0.intrcpt, 'lower'));
           x0.arch_param;
           x0.garch_param;
           vech(chol(x0.war_matrix, 'lower'));
           x0.df_F_1;
           x0.df_F_2
    ];           
end
A = [zeros(1,k*(k+1)/2) 1 1 zeros(1,k*(k+1)/2+2)]; % arch_param + garch_param < 1
b = [1];
Aeq = [];
beq = [];
lb = [ -inf(k*(k+1)/2,1);
       0;
       0;
       -inf(k*(k+1)/2,1);
       k;
       k
];
ub = [ inf(k*(k+1)/2,1);
       1;%inf(p,1);
       1;%inf(q,1);
       inf(k*(k+1)/2,1);
       inf;
       inf
];
nonlcon = [];
options = optimoptions('fmincon',...
    'Display', 'iter-detailed',... % 'iter-detailed',... % 'PlotFcn', @optimplotx,... 
    'SubproblemAlgorithm', 'cg',... % This works better :) 'DiffMinChange', 1e-8,...
    'MaxFunctionEvaluations', 1e5,...
    'UseParallel', false...
);

warning('off','MATLAB:nearlySingularMatrix') %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[eparam,~,exitflag,optimoutput,lambda,grad,hessian] = fmincon(...
    fun,x0,A,b,Aeq,beq,lb,ub,nonlcon,options);
warning('on','MATLAB:nearlySingularMatrix') %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Output
[intrcpt, arch_param, garch_param, war_matrix, df_F_1, df_F_2] = ...
    param_trans(eparam, k);

[ nLogL, logLcontr, Sigma_, Omega_ ] = ...
    cancF_likeRec(intrcpt, arch_param, garch_param, ...
                  war_matrix, df_F_1, df_F_2, data_RC);

eparam = struct(...
    'intrcpt', intrcpt,...
    'arch_param', arch_param,...
    'garch_param', garch_param,...
    'war_matrix', war_matrix,...
    'df_F_1', df_F_1,...
    'df_F_2', df_F_2...
);
logL = struct('logL', -nLogL, 'logLcontr', logLcontr);
fit = struct('Sigma_', Sigma_(:,:,1:T), 'Omega_', Omega_(:,:,1:T));
fcst = struct('Sigma_', Sigma_(:,:,T+1:end), 'Omega_', Omega_(:,:,T+1:end));

end

%% Wrapper function for likeRec to take parameters input in single vector
% Needed for optimization.
function [nLogL,logLcontr,Sigma_,Omega_] = likeRec_wrapper(param, data_RC) 

[k,~,~] = size(data_RC);
if numel(param) ~= (k*(k+1) + 4)
    error('Number of elements in param does not match rest of inputs.')
end

[intrcpt, arch_param, garch_param, war_matrix, df_F_1, df_F_2] = param_trans(param, k);

[ nLogL, logLcontr, Sigma_, Omega_ ] = cancF_likeRec( ...
    intrcpt, arch_param, garch_param, war_matrix, df_F_1, df_F_2, data_RC );
end

%% Function to transform parameters from vector form to explicit form
function [intrcpt, arch_param, garch_param, war_matrix, df_F_1, df_F_2] = param_trans(param, k)

intrcpt = ivech(param(1 : k*(k+1)/2), 'lower');
intrcpt = intrcpt*intrcpt';

arch_param = param(k*(k+1)/2+1);
garch_param = param(k*(k+1)/2+2);

war_matrix = ivech(param(k*(k+1)/2+3 : k*(k+1)/2+2+k*(k+1)/2), 'lower');
war_matrix = war_matrix*war_matrix';

df_F_1 = param(k*(k+1) + 3);
df_F_2 = param(k*(k+1) + 4);
end
