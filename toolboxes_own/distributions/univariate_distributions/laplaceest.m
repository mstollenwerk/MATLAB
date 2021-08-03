function [eparam, tstats, logL] = laplaceest(dta)
%LAPLACEEST MLE of univariate Laplace distribution
%
% REFERENCES:
%      [1] Kotz, Kozubowski and Podgorski (2001), p. 71,74, 6.5.3

narginchk(1,1);
if size(dta,1) < size(dta,2)
    dta = dta';
end
N = size(dta,1);
%% MLE
eparam.theta_ = median(dta);
eparam.s = mean( abs( dta - eparam.theta_ ) );
%% tstats
tstats.theta_ = eparam.theta_/eparam.s^2;
tstats.s = eparam.s;
%% nLogL, logLcontr and eparam
logLcontr = log(laplacepdf(eparam.theta_, eparam.s, dta));
nLogL = -sum(logLcontr);

aic = 2*nLogL + 4;
bic = 2*nLogL + log(N)*2;
logL = struct(...
    'logL', -nLogL, ...
    'logLcontr', logLcontr, ...
    'bic', bic, ...
    'aic', aic ...
);
end

