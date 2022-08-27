clear
clc
load('C:\Users\Stollenwerk\OneDrive - bwedu\Research\Projekte\Laufend\rc_distributions_matter\data\rc.mat')
rv_used = squeeze(rc(1,1,:));
m=63; % 21, 63, 126, 252

x0 = [0.007;0.85; 0.01; 0.07; 0.91 ]; % m = 63
options = optimoptions(@fminunc, "Algorithm", "active-set", "Display", "iter-detailed", 'MaxFunEval', 1e4, 'MaxIter', 1e4);
eparam = fminunc('likelihood_mf2_garch', x0, options, rv_used, m);

[ e, h, tau, V_m ] = mf2_garch_core(eparam, rv_used, m);

