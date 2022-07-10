function R = matvrnd(dist, Omega_, dfs, T)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
p = size(Omega_,1);

if strcmp( dist, 'Wish' )
    if length(dfs) ~= 1
        error('length(dfs) wrong.')
    end    
    n = dfs(1);
    R = matvWishrnd(Omega_,n,T);
elseif strcmp( dist, 'iWish' )
    if length(dfs) ~= 1
        error('length(dfs) wrong.')
    end    
    nu = dfs(1);
    R =  matviWishrnd(Omega_,nu,T);   
elseif strcmp( dist, 'tWish' )
    if length(dfs) ~= 2
        error('length(dfs) wrong.')
    end    
    n = dfs(1); 
    nu = dfs(2); 
    R =  matvtWishrnd(Omega_,n,nu,T);    
elseif strcmp( dist, 'itWish' )
    if length(dfs) ~= 2
        error('length(dfs) wrong.')
    end        
    n = dfs(1); 
    nu = dfs(2);
    R =  matvitWishrnd(Omega_,n,nu,T);
elseif strcmp( dist, 'F' )
    if length(dfs) ~= 2
        error('length(dfs) wrong.')
    end    
    n = dfs(1); 
    nu = dfs(2);
    R =  matvFrnd(Omega_,n,nu,T);    
elseif strcmp( dist, 'Riesz' )
    if length(dfs) ~= p
        error('length(dfs) wrong.')
    end        
    n = dfs(1 : p);
    R = matvRieszrnd(Omega_,n,T);
elseif strcmp( dist, 'iRiesz2' )
    if length(dfs) ~= p
        error('length(dfs) wrong.')
    end        
    nu = dfs(1 : p); 
    R = matviRiesz2rnd(Omega_,nu,T); 
elseif strcmp( dist, 'tRiesz' )
    if length(dfs) ~= p+1
        error('length(dfs) wrong.')
    end        
    n = dfs(1 : p); 
    nu = dfs(p + 1);
    R = matvtRieszrnd(Omega_,n,nu,T);
elseif strcmp( dist, 'itRiesz2' )
    if length(dfs) ~= p+1
        error('length(dfs) wrong.')
    end    
    n = dfs(1); 
    nu = dfs(2 : p + 1);
    R =  matvitRiesz2rnd(Omega_,n,nu,T); 
elseif strcmp( dist, 'FRiesz' )
    if length(dfs) ~= p+p
        error('length(dfs) wrong.')
    end    
    n = dfs(1 : p); 
    nu = dfs(p + 1 : p + p);
    R =  matvFRieszrnd(Omega_,n,nu,T); 
elseif strcmp( dist, 'iFRiesz2' )
    if length(dfs) ~= p+p
        error('length(dfs) wrong.')
    end    
    n = dfs(1 : p); 
    nu = dfs(p + 1 : p + p);
    R =  matviFRiesz2rnd(Omega_,n,nu,T); 
end


end