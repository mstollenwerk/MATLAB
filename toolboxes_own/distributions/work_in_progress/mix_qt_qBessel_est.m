function [ eparam, tstats, logL, optimoutput ] = mix_qt_qBessel_est( data_mat, varargin )
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
fun = @( param ) like_wrapper(param, data_mat);
% x0 ----------------------------------------------------------------------
if isempty(x0) || isstruct(x0) 
    if isfield(x0,'b_df_1')
        b_df_1 = x0.b_df_1;
    else
        b_df_1 = 2*p;
    end
    if isfield(x0,'b_df_2')
        b_df_2 = x0.b_df_2;
    else
        b_df_2 = 2*p;
    end 
    if isfield(x0,'b_df_3')
        b_df_3 = x0.b_df_3;
    else
        b_df_3 = 2*p;
    end  
    if isfield(x0,'b_Sigma')
        b_Sigma = x0.b_Sigma;
    else
        b_Sigma = eye(p); %mean(data_mat,3);
    end     
    if isfield(x0,'t_df_1')
        t_df_1 = x0.t_df_1;
    else
        t_df_1 = 2*p;
    end
    if isfield(x0,'t_df_2')
        t_df_2 = x0.t_df_2;
    else
        t_df_2 = 2*p;
    end     
    if isfield(x0,'b_Sigma')
        t_Sigma = x0.t_Sigma;
    else
        t_Sigma = eye(p); %mean(data_mat,3);
    end
    if isfield(x0,'mixture')
        mixture = x0.mixture;
    else
        mixture = .5;
    end         
    x0 = [ vech(chol(b_Sigma, 'lower')); vech(chol(t_Sigma, 'lower')); b_df_1; b_df_2; b_df_3; t_df_1; t_df_2; mixture]';
    try
        fun(x0);
    catch
        error('x0 is not valid.')
    end
else
    error('x0 must be empty or struct.')
end
% Optimoptions-------------------------------------------------------------
options = optimoptions('fmincon',...
    'MaxFunctionEvaluations', 1e5 ...
);
options = optimoptions(options, varargin{2:end});
% Restrictions-------------------------------------------------------------
A = [ zeros(1,2*p_) -p/2 -1 0 0 0 0 ]; %
b = [-eps];
Aeq = [];
beq = [];
lb = [ -inf(2*p_,1);
       p;
       -p*N/2;
       0;
       4;
       1;
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
[nLogL, logLcontr, b_Sigma, t_Sigma, b_df_1, b_df_2, b_df_3, t_df_1, t_df_2, mixture ] = like_wrapper(eparam, data_mat);

eparam = struct( ...
    'b_Sigma', b_Sigma, ...
    't_Sigma', t_Sigma, ...
    'b_df_1', b_df_1, ...
    'b_df_2', b_df_2, ...    
    'b_df_3', b_df_3, ...   
    't_df_1', t_df_1, ...
    't_df_2', t_df_2, ... 
    'mixture', mixture, ...
    'all', eparam ...
);

aic = 2*nLogL + 2*numel(x0);
bic = 2*nLogL + log(N)*numel(x0);
logL = struct(...
    'logL', -nLogL, ...
    'logLcontr', logLcontr, ...
    'bic', bic, ...
    'aic', aic ...
);
tstats = struct(...
    'b_Sigma_', tstats(1 : p_),...
    't_Sigma_', tstats(p_ + 1 : 2*p_),...
    'b_df_1', tstats(2*p_ + 1), ...
    'b_df_2', tstats(2*p_ + 2), ...    
    'b_df_3', tstats(2*p_ + 3), ...        
    't_df_1', tstats(2*p_ + 4), ...     
    't_df_2', tstats(2*p_ + 5), ...    
    'mixture', tstats(2*p_ + 6), ...    
    'all', tstats ...
);
%% Likelihood wrapper function
function [nLogL, logLcontr, b_Sigma, t_Sigma, b_df_1, b_df_2, b_df_3, t_df_1, t_df_2, mixture ] = like_wrapper(param, data_mat) 

b_Sigma = ivechchol(param(1:p_));
t_Sigma = ivechchol(param(p_+1:2*p_));
b_df_1 = param(2*p_+1);
b_df_2 = param(2*p_+2);
b_df_3 = param(2*p_+3);
t_df_1 = param(2*p_+4);
t_df_2 = param(2*p_+5);
mixture = param(2*p_+6);

[ ~, b_logLcontr ] = qBessellike( ...
    data_mat, ...
    ones(1,N)*b_df_1, ...
    ones(1,N)*b_df_2, ...
    ones(1,N)*b_df_3, ...
    repmat(b_Sigma,1,1,N) ...
);

[ ~, t_logLcontr ] = qtlike( ...
    t_Sigma, ...
    t_df_1, ...
    t_df_2, ...    
    data_mat ...
);

logLcontr = log(mixture*exp(b_logLcontr) + (1-mixture)*exp(t_logLcontr));
nLogL = -sum(logLcontr);

end
end
