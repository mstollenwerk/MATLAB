function clusterdummy()
load('everythingUneed') 

parpool('local', str2num(getenv('SLURM_CPUS_PER_TASK')) ) % SLURM Cluster ,'AttachedFiles',{'everythingUneed1.mat',{everythingUneed2.m'}
parpool('local',str2num(getenv('MOAB_PROCCOUNT'))); % MOAB Cluster

poolobj = gcp('nocreate'); % If no pool, do not create new one.
if isempty(poolobj)
       poolSize = 0;
       %error('parallel:demo:poolClosed', ...
       %'This demo needs an open MATLAB pool to run.');
else
       poolSize = poolobj.NumWorkers;
end

%
% Code Here
%

save('everythingUneed','everythingUneed')
delete(gcp);
end