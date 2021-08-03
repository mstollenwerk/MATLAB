function [eparam, tstats, logL, optimoutput] = matvFGammamixest(X, x0, varargin)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

%% Input Checking
% Will be added later
[p,~,N] = size(X);
p_ = p*(p+1)/2;
%% Optimization
if isempty(x0)
    x0 = [vech(chol(eye(p),'lower'));2*p;2*p;5];
end
% Optimoptions-------------------------------------------------------------
options = optimoptions('fmincon',...
    'OutputFcn', @outfun_rec_x, ...        
    'Display', 'iter-detailed' ...
);
options = optimoptions(options, varargin{:});
% Restrictions-------------------------------------------------------------
lb = [ -inf(p_,1); p-1; 0; 0];
% Optimization-------------------------------------------------------------
tic

% Set up shared variables with OUTFUN
history.x = [];
history.fval = [];

warning('off')
[eparam,~,exitflag,optimoutput,grad,hessian] = ...
    fmincon( ...
        @(all_param) matvFGammamixlike([],[],[],[],X,all_param), ...
        x0,[],[],[],[],lb,[],[],options ...
    );
warning('on')

optimoutput.estimation_time = toc;
optimoutput.exitflag = exitflag;
optimoutput.grad = grad;
optimoutput.hessian = hessian;

optimoutput.history.x = history.x;
optimoutput.history.fval = history.fval;

[ nLogL, logLcontr, eparam ] = matvFGammamixlike([],[],[],[],X,eparam);
logL.logL = -nLogL;
logL.logLcontr = logLcontr;
logL.aic = 2*nLogL + 2*numel(x0);
logL.bic = 2*nLogL + log(N)*numel(x0);

% Tstats-------------------------------------------------------------------
%[VCV,A,B,scores,hess,gross_scores] = robustvcv(fun, eparam, 3);
% [VCV,scores,gross_scores] = ...
%     vcv(...
%         @(all_param) matvFGammamixlike([],[],[],[],X,all_param), ...
%         eparam.all ...
%     );
% tstats = eparam.all./sqrt(diag(VCV)');
tstats = NaN;

%% Outfunction to record values during optimization
    % It has to be a local function.
    function stop = outfun_rec_x(x,optimValues,state)
    % This is an output function for optimization algos. It records the points,
    % fvals and directions. 
    %
    % The history.x, history.fval and searchdir variables have to be set up in
    % the main code prior to optimization.

    stop = false;

    switch state
       case 'init'
           history.fval = [history.fval; optimValues.fval];
           history.x = [history.x; x];
       case 'iter'
           history.fval = [history.fval; optimValues.fval];
           history.x = [history.x; x];
       otherwise
    end
    end
end

