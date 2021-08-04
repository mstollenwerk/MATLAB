function [eparam, logL, optimoutput] = matvLaplaceFest(X, x0, varargin)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

%% Input Checking
% Will be added later
[k,~,~] = size(X);
%% Optimization
if isempty(x0)
    x0 = [vech(chol(eye(k),'lower'));2*k;2*k;1];
end
% Optimoptions-------------------------------------------------------------
options = optimoptions('fmincon',...
    'SubproblemAlgorithm', 'cg',...
    'StepTolerance', 1e-05, ...    
    'OutputFcn', @outfun_rec_x, ...        
    'MaxFunctionEvaluations', 1e5, ...
    'Display', 'iter-detailed' ...
);
options = optimoptions(options, varargin{:});
% Restrictions-------------------------------------------------------------
lb = [ -inf(k*(k+1)/2,1); k; k; 0];
% Optimization-------------------------------------------------------------
tic

% Set up shared variables with OUTFUN
history.x = [];
history.fval = [];


warning('off')
[eparam,~,exitflag,optimoutput,grad,hessian] = ...
    fmincon( ...
        @(all_param) matvLaplaceFlike([],[],[],[],X,all_param), ...
        x0,[],[],[],[],lb,[],[],options ...
    );
warning('on')

optimoutput.estimation_time = toc;
optimoutput.exitflag = exitflag;
optimoutput.grad = grad;
optimoutput.hessian = hessian;

optimoutput.history.x = history.x;
optimoutput.history.fval = history.fval;

[ nLogL, logLcontr, eparam ] = matvLaplaceFlike([],[],[],[],X,eparam);
logL.logL = -nLogL;
logL.logLcontr = logLcontr;

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

