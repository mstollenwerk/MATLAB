function [eparam,optimoutput] = matvRieszEstGradHess(R,varargin)

p = size(R,1);

obj_fun = @(n) obj_fun_wrapper(mean(R,3), n, R);

x0 = 2*p*ones(p,1);

options = optimoptions('fminunc','Algorithm','trust-region',...
    'SpecifyObjectiveGradient',true,'HessianFcn','objective','UseParallel',true,'Display','iter-detailed');

tic
warning('off')
[eparam,fval,exitflag,optimoutput] = fminunc(obj_fun,x0,options);
warning('on')
optimoutput.exitflag = exitflag;
optimoutput.fval = fval;
optimoutput.estimation_time = toc;

end

function [nLogL, g, H] = obj_fun_wrapper(Sigma_, n, R)

Omega_ = matvStandardize('Riesz', Sigma_, n);
try
    [ nLogL, g, H ] = matvRieszlikeGradHess( Omega_, n, R );
catch
    nLogL = inf;
    g = NaN;
    H = NaN;
end

end