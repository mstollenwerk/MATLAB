function [X,optimoutput] = my_fminunc(FUN,X,varargin)
%MY_FMINUNC Summary of this function goes here
%   Detailed explanation goes here

options = optimoptions('fminunc',... 
    'OutputFcn', @outfun_rec_x, ...        
    'Display', 'iter-detailed' ...
);
options = optimoptions(options, varargin{:});

tic

% Set up shared variables with OUTFUN
history.x = [];
history.fval = [];

warning('off')
[X,~,exitflag,optimoutput,grad,hessian] = fminunc(FUN,X,options);
warning('on')

% % If solver stops at larger value than the smallest value evalutated during
% % optimization, let it start from that smallest value again.
% while history.fval(end) ~= min(history.fval)
%     X = history.x(history.fval==min(history.fval),:);
% 
%     history.x = [];
%     history.fval = [];
% 
%     [X,~,exitflag,optimoutput,grad,hessian] = fminunc(FUN,X,options);
% end

optimoutput.estimation_time = toc;
optimoutput.exitflag = exitflag;
optimoutput.grad = grad;
optimoutput.hessian = hessian;

optimoutput.history.x = history.x;
optimoutput.history.fval = history.fval;

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

