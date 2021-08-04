function [ eparam, tstats, logL, optimoutput ] = qBesselest( data_mat, varargin )
%BESSELWISHEST
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
%% Error checking
narginchk(1,inf);
if nargin==1
    x0 = [];
else
    x0 = varargin{1};
end
[p,p2,N] = size(data_mat);
if p~=p2
    error('Covariance matrix data must be square.')
end
p_ = p*(p+1)/2;
%% Optimization
% Objective Function ------------------------------------------------------
fun = @( param ) like_wrapper(param, data_mat, N);
% x0 ----------------------------------------------------------------------
if isempty(x0) || isstruct(x0) 
    if isfield(x0,'df_1')
        df_1 = x0.df_1;
    else
        df_1 = 2*p;
    end
    if isfield(x0,'df_2')
        df_2 = x0.df_2;
    else
        df_2 = 2*p;
    end 
    if isfield(x0,'df_3')
        df_3 = x0.df_3;
    else
        df_3 = 2*p;
    end      
    if isfield(x0,'Sigma_')
        Sigma_ = x0.Sigma_;
    else
        Sigma_ = eye(p); %mean(data_mat,3);
    end      
    x0 = [ vech(chol(Sigma_, 'lower')); df_1; df_2; df_3 ]'; % k*(k+1)/2 + 3
    try
        fun(x0);
    catch
        error('x0 is not valid.')
    end
else
    error('x0 must be empty or struct.')
end
% Restrictions-------------------------------------------------------------
% Optimoptions-------------------------------------------------------------
options = optimoptions('fmincon',...
    'MaxFunctionEvaluations', 1e5 ...
);
options = optimoptions(options, varargin{2:end});
% Restrictions-------------------------------------------------------------
A = [ zeros(1,p_) -p/2 -1 0 ]; %
b = [-eps];
Aeq = [];
beq = [];
lb = [ -inf(p_,1);
       p;
       -p*N/2;
       0
];
ub = [];
nonlcon = [];
% Optimization-------------------------------------------------------------
warning('off','MATLAB:nearlySingularMatrix') %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic
[eparam,~,exitflag,optimoutput,grad,hessian] = fmincon(...
    fun,x0,A,b,Aeq,beq,lb,ub,nonlcon,options);
optimoutput.estimation_time = toc;
optimoutput.exitflag = exitflag;
optimoutput.grad = grad;
optimoutput.hessian = hessian;
warning('on','MATLAB:nearlySingularMatrix') %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Tstats-------------------------------------------------------------------
%[VCV,A,B,scores,hess,gross_scores] = robustvcv(fun, eparam, 3);
[VCV,scores,gross_scores] = vcv(fun, eparam);
tstats = eparam./sqrt(diag(VCV)');
% Output Creation----------------------------------------------------------
[nLogL, logLcontr, Sigma_, df_1, df_2, df_3 ] = like_wrapper(eparam, data_mat, N);

eparam = struct( ...
    'Sigma_', Sigma_, ...
    'df_1', df_1, ...
    'df_2', df_2, ...    
    'df_3', df_3, ...       
    'all', eparam ...
);

aic = 2*nLogL + 2*(numel(x0));
bic = 2*nLogL + log(N)*(numel(x0));
logL = struct(...
    'logL', -nLogL, ...
    'logLcontr', logLcontr, ...
    'bic', bic, ...
    'aic', aic ...
);
tstats = struct(...
    'Sigma_', tstats(1 : p_),...
    'df_1', tstats(p_ + 1), ...
    'df_2', tstats(p_ + 2), ...    
    'df_3', tstats(p_ + 3) ...        
);
%% Likelihood wrapper function
function [nLogL, logLcontr, Sigma_, df_1, df_2, df_3 ] = like_wrapper(param, data_mat, N) 

Sigma_ = ivech(param(1:p_), 'lower');
Sigma_ = Sigma_*Sigma_';
df_1 = param(p_+1);
df_2 = param(p_+2);
df_3 = param(p_+3);

[ nLogL, logLcontr ] = qBessellike( ...
    data_mat, ...
    ones(1,N)*df_1, ones(1,N)*df_2, ones(1,N)*df_3, ...
    repmat(Sigma_,1,1,N) ...
);

end
end
