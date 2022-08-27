function [eparam, logL , fit] = mf2heavyGarchEstim(r,RC,varargin)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
obj_fun = @(param) mf2heavyGarchLikeRec(param,r,RC);

options = optimoptions('fminunc', 'Display', 'iter-detailed', 'MaxFunEval', 1e4);

if nargin ==3
    x0 = varargin{1};
else
    x0 = [diag(mean(RC,3))', .05, .9, .05, .9, .05, .9]';
end

eparam = fminunc(obj_fun,x0,options);

[nLogL,logLcontr,h,z,Corr,H] = obj_fun(eparam);

logL.nLogL = nLogL;
logL.logLcontr = logLcontr;

fit.h = h;
fit.z = z;
fit.Corr = Corr;
fit.H = H;

end

