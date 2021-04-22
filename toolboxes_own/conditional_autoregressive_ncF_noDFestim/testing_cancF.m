clear
clc

p = 2;
mcrep = 100;
T = 2000;

param.intrcpt = randi(15,2);
param.intrcpt = param.intrcpt*param.intrcpt'/100;
param.arch_param = .075;
param.garch_param = .9;
param.war_matrix = randi(10,2);
param.war_matrix = param.war_matrix*param.war_matrix'/100;

parpool(25)
parfor ii = 1:mcrep
    [ data_RC{ii}, Sigma_{ii}, Omega_{ii}] = toymodel_rnd( param, T );
    [eparam{ii}, logL{ii}, fit_fcst{ii}, optimoutput{ii}] = toymodel_estim(data_RC{ii}, []);
end

eparam_ = NaN(mcrep,p*(p+1)+2);
eintrcpt = NaN(p,p,mcrep);
ewar_matrix = NaN(p,p,mcrep);
for ii = 1:mcrep
    eintrcpt(:,:,ii) = eparam{ii}.intrcpt;
    ewar_matrix(:,:,ii) = eparam{ii}.war_matrix;
    eparam_(ii,1:3) = vech(chol(eparam{ii}.intrcpt, 'lower'));
    eparam_(ii,4) = eparam{ii}.arch_param;
    eparam_(ii,5) = eparam{ii}.garch_param;
    eparam_(ii,6:8) = vech(chol(eparam{ii}.war_matrix, 'lower'));
end
