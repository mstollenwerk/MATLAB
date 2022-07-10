clear
clc

if isempty(gcp('nocreate'))
    parpool(maxNumCompThreads,'IdleTimeout', 300)
end

p = 3;
N = 1e4;

dists = {'Wish'
        'iWish'
        'tWish'
        'itWish'
        'F'
        'Riesz'
        'iRiesz2'
        'tRiesz'
        'itRiesz2'
        'FRiesz'
        'iFRiesz2'};
    
lengthdf= [1, 1, 2, 2, 2, p, p, p+1, p+1, 2*p, 2*p];

Sig = rndpd(p);
n = chi2rnd(p+10,p,1);
nu = chi2rnd(p+10,p,1);
dfs{1} = n(1);
dfs{2} = nu(1);
dfs{3} = [n(1); nu(1)];
dfs{4} = [n(1); nu(1)];
dfs{5} = [n(1); nu(1)];
dfs{6} = n;
dfs{7} = nu;
dfs{8} = [n;nu(1)];
dfs{9} = [n(1);nu];
dfs{10} = [n;nu];
dfs{11} = [n;nu];


parfor dist_id = 1:length(dists)
    dist = dists(dist_id);       
    dfs_dist = dfs{dist_id};
    
    R = NaN(p,p,N);
    Om = matvStandardize(dist,Sig,dfs_dist);
    R(:,:,1:N/2) = matvrnd(dist,Om,dfs_dist,N/2);
    R(:,:,N/2+1:N) = matvsrnd(dist,Sig,dfs_dist,N/2);

    [ eparam{dist_id}, tstats{dist_id}, logL{dist_id}, optimoutput{dist_id}] = ...
        matvest(R,dist,[],'noperm','UseParallel',true);
    
end