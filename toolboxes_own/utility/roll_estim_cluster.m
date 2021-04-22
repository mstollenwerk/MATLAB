function [eparam, eparam_tstats, logL, fit_fcst, optimoutput] = roll_estim_cluster(Tstart, Tend, wndw)
%ROLL_ESTIM Estimates a model on rolling window of data from Tstart to Tend
%   Detailed explanation goes here

%% Input Checking
if Tstart - wndw < 0 
    error('Tstart is smaller than window')
end
%% Get Data, determine T and inputdata to estimation
load('C:\Users\Micha\OneDrive - bwstaff\Research\minute_data\rc_data\rcov_aaba_ba.mat')
data_y = ret_otc_aaba_ba;
data_RC = rc_aaba_ba;
T = size(data_RC,3);

%% Rolling estimation| Change model here
eparam = cell(T,1);
eparam_tstats = cell(T,1);
logL = cell(T,1);
fit_fcst = cell(T,1);
optimoutput = cell(T,1);

pc = parcluster('local');
pc.JobStorageLocation = getenv('TMPDIR');
num_workers = str2num(getenv('MOAB_PROCCOUNT'));
parpool(pc, num_workers);
pc.NumWorkers=num_workers;

parfor tt = Tstart:Tend
    [eparam{tt}, eparam_tstats{tt}, logL{tt}, fit_fcst{tt}, optimoutput{tt}] = ...
        heavy_gas_t_F_estim(data_y(tt-wndw+1:tt,:), data_RC(:,:,tt-wndw+1:tt), 1, 1, [], [] );
end

save([num2str(Tstart),'_',num2str(Tend),'.mat'], 'eparam', 'eparam_tstats', 'logL', 'fit_fcst', 'optimoutput', 'dates', 'ret_otc_aaba_ba', 'rc_aaba_ba')
end

