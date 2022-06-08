function R = matvsrnd(dist, Sigma_, dfs, T)
%MATVSRND Random matrices following a set of standardized distributions.
%   Detailed explanation goes here
p = size(Sigma_,1);

if strcmp( dist, 'Wish' )
    if length(dfs) ~= 1
        error('length(dfs) wrong.')
    end    
    n = dfs(1);
    R = matvWishrnd(Sigma_,n,T)./n;
elseif strcmp( dist, 'iWish' )
    if length(dfs) ~= 1
        error('length(dfs) wrong.')
    end    
    nu = dfs(1);
    R =  matviWishrnd(Sigma_,nu,T).*(nu-p-1);   
elseif strcmp( dist, 'tWish' )
    if length(dfs) ~= 2
        error('length(dfs) wrong.')
    end    
    n = dfs(1); 
    nu = dfs(2); 
    R =  matvtWishrnd(Sigma_,n,nu,T).*((nu-2)/n/nu);    
elseif strcmp( dist, 'itWish' )
    if length(dfs) ~= 2
        error('length(dfs) wrong.')
    end        
    n = dfs(1); 
    nu = dfs(2);
    R =  matvitWishrnd(Sigma_,n,nu,T).*(n-p-1);
elseif strcmp( dist, 'F' )
    if length(dfs) ~= 2
        error('length(dfs) wrong.')
    end    
    n = dfs(1); 
    nu = dfs(2);
    R =  matvFrnd(Sigma_,n,nu,T).*((nu-p-1)/n);    
elseif strcmp( dist, 'Riesz' )
    if length(dfs) ~= p
        error('length(dfs) wrong.')
    end        
    n = dfs(1 : p);
    R =  matvsRieszrnd(Sigma_,n,T);
elseif strcmp( dist, 'iRiesz2' )
    if length(dfs) ~= p
        error('length(dfs) wrong.')
    end        
    nu = dfs(1 : p); 
    R =  matvsiRiesz2rnd(Sigma_,nu,T); 
elseif strcmp( dist, 'tRiesz' )
    if length(dfs) ~= p+1
        error('length(dfs) wrong.')
    end        
    n = dfs(1 : p); 
    nu = dfs(p + 1);
    R =  matvstRieszrnd(Sigma_,n,nu,T);
elseif strcmp( dist, 'itRiesz2' )
    if length(dfs) ~= p+1
        error('length(dfs) wrong.')
    end    
    n = dfs(1); 
    nu = dfs(2 : p + 1);
    R =  matvsitRiesz2rnd(Sigma_,n,nu,T); 
elseif strcmp( dist, 'FRiesz' )
    if length(dfs) ~= p+p
        error('length(dfs) wrong.')
    end    
    n = dfs(1 : p); 
    nu = dfs(p + 1 : p + p);
    R =  matvsFRieszrnd(Sigma_,n,nu,T); 
elseif strcmp( dist, 'iFRiesz2' )
    if length(dfs) ~= p+p
        error('length(dfs) wrong.')
    end    
    n = dfs(1 : p); 
    nu = dfs(p + 1 : p + p);
    R =  matvsiFRiesz2rnd(Sigma_,n,nu,T); 
end


end