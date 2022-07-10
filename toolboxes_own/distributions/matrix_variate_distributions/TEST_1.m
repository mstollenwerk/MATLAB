clear
clc

p = 3;
N = 100;

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

for dist_id = 1:length(dists)
    dist = dists(dist_id);
    
    Om = rndpd(p);
    dfs = chi2rnd(p+10,lengthdf(dist_id),1);

    R = NaN(p,p,N);
    Sig = matvEV(dist,Om,dfs);
    R(:,:,1:N/2) = matvrnd(dist,Om,dfs,N/2);
    R(:,:,N/2+1:N) = matvsrnd(dist,Sig,dfs,N/2);

    [~,~,score] = matvLogLike(dist,Om,dfs,R);
    [~,~,score_s] = matvsLogLike(dist,Sig,dfs,R);
    
    if isfield(score,'n')
        
        figure
        plot(score.n,'LineStyle','none','Marker','o')
        hold on
        plot(score_s.n_originalpdf,'LineStyle','none','Marker','o')
        plot(repmat(mean(score.n),N,1))
        title(strcat(dists(dist_id),', n'))
        
        figure
        plot(score.n_scaled,'LineStyle','none','Marker','o')
        hold on
        plot(score_s.n_originalpdf_scaled,'LineStyle','none','Marker','o')
        plot(repmat(mean(score.n_scaled),N,1))
        title(strcat(dists(dist_id),', n, scaled'))
        
    end
    
    if isfield(score,'nu')
        
        figure
        plot(score.nu,'LineStyle','none','Marker','o')
        hold on
        plot(score_s.nu_originalpdf,'LineStyle','none','Marker','o')
        plot(repmat(mean(score.nu),N,1))
        title(strcat(dists(dist_id),', nu'))
        
        figure
        plot(score.nu_scaled,'LineStyle','none','Marker','o')
        hold on
        plot(score_s.nu_originalpdf_scaled,'LineStyle','none','Marker','o')
        plot(repmat(mean(score.nu_scaled),N,1))
        title(strcat(dists(dist_id),', nu, scaled'))
        
    end    
end