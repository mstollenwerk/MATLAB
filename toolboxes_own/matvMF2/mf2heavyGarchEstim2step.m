function [eparam, logL , fit] = mf2heavyGarchEstim2step(r,RC,varargin)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

[p,~,T] = size(RC);

options = optimoptions('fminunc', 'Display', 'iter-detailed', 'MaxFunEval', 1e4);

%% Step 1
mStep1 = 6;

obj_funStep1 = @(param) mf2heavyGarchLikeRecStep1(param,r,RC,mStep1);

if nargin ==3
    x0Step1 = varargin{1}(1:end-2);
else
    x0Step1 = [0.05*ones(1,p), 0.2, 0.6, 0.2, 0.6]';
end

eparamStep1 = fminunc(obj_funStep1,x0Step1,options);
[nLogLStep1,logLcontrStep1,h,tau] = obj_funStep1(eparamStep1);
sig = h.*tau;

%% Step 2
mStep2 = 1;

rStan = r./sqrt(sig);

RCorr = NaN(p,p,T);
for tt = 1:T
    rsdevs = sqrt(diag(RC(:,:,tt)));
    RCorr(:,:,tt) = RC(:,:,tt)./(rsdevs*rsdevs');
end

obj_funStep2 = @(param) mf2heavyGarchLikeRecStep2(param,rStan,RCorr,mStep2);

if nargin ==3
    x0Step2 = varargin{1}(end-1:end);
else
    x0Step2 = [.05, .9]';
end

eparamStep2 = fminunc(obj_funStep2,x0Step1,options);

[nLogLStep2,logLcontrStep2,RCorrE] = obj_funStep2(eparamStep2);

%% Output
eparam.all = [eparamStep1, eparamStep2];
eparam.intrcptTau = eparamStep1(1:p);
eparam.aTau = eparamStep1(p+1);
eparam.bTau = eparamStep1(p+2);
eparam.aS = eparamStep1(p+3);
eparam.bS = eparamStep1(p+4);
eparam.aCorr = eparamStep2(1);
eparam.bCorr = eparamStep2(2);

logL.nLogL = nLogLStep1 + nLogLStep2;
logL.logLcontr = logLcontrStep1 + logLcontrStep2;

fit.h = h;
fit.tau = tau;
fit.sig = sig;
fit.Corr = RCorrE;

end

