function A = matvrnd(dist, Omega_, dfs, T)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
k = size(Omega_,1);

if strcmp( dist, 'Wish' )
    if length(dfs) ~= 1
        error('length(dfs) wrong.')
    end    
    n = dfs(1);
    A =  matvWishrnd(Omega_,n,T);
elseif strcmp( dist, 'iWish' )
    if length(dfs) ~= 1
        error('length(dfs) wrong.')
    end    
    n = dfs(1);
    A =  matviWishrnd(Omega_,n,T);   
elseif strcmp( dist, 'tWish' )
    if length(dfs) ~= 2
        error('length(dfs) wrong.')
    end    
    n = dfs(1); 
    nu = dfs(2); 
    A =  matvtWishrnd(Omega_,n,nu,T);    
elseif strcmp( dist, 'itWish' )
    if length(dfs) ~= 2
        error('length(dfs) wrong.')
    end        
    n = dfs(1); 
    nu = dfs(2);
    A =  matvitWishrnd(Omega_,n,nu,T);
elseif strcmp( dist, 'F' )
    if length(dfs) ~= 2
        error('length(dfs) wrong.')
    end    
    n = dfs(1); 
    nu = dfs(2);
    A =  matvFrnd(Omega_,n,nu,T);    
elseif strcmp( dist, 'Riesz' )
    if length(dfs) ~= k
        error('length(dfs) wrong.')
    end        
    n = dfs(1 : k);
    A =  matvRieszrnd(Omega_,n,T);    
elseif strcmp( dist, 'Riesz2' )
    if length(dfs) ~= k
        error('length(dfs) wrong.')
    end        
    n = dfs(1 : k);
    A =  matvRiesz2rnd(Omega_,n,T);        
elseif strcmp( dist, 'iRiesz' )
    if length(dfs) ~= k
        error('length(dfs) wrong.')
    end        
    n = dfs(1 : k);
    A =  matviRieszrnd(Omega_,n,T); 
elseif strcmp( dist, 'iRiesz2' )
    if length(dfs) ~= k
        error('length(dfs) wrong.')
    end        
    n = dfs(1 : k); 
    A =  matviRiesz2rnd(Omega_,n,T); 
elseif strcmp( dist, 'tRiesz' )
    if length(dfs) ~= k+1
        error('length(dfs) wrong.')
    end        
    n = dfs(1 : k); 
    nu = dfs(k + 1);
    A =  matvtRieszrnd(Omega_,n,nu,T); 
elseif strcmp( dist, 'tRiesz2' )
    if length(dfs) ~= k+1
        error('length(dfs) wrong.')
    end    
    n = dfs(1 : k); 
    nu = dfs(k + 1); 
    A =  matvtRiesz2rnd(Omega_,n,nu,T); 
elseif strcmp( dist, 'itRiesz' )
    if length(dfs) ~= k+1
        error('length(dfs) wrong.')
    end    
    n = dfs(1 : k); 
    nu = dfs(k + 1);
    A =  matvitRieszrnd(Omega_,n,nu,T); 
elseif strcmp( dist, 'itRiesz2' )
    if length(dfs) ~= k+1
        error('length(dfs) wrong.')
    end    
    n = dfs(1 : k); 
    nu = dfs(k + 1);
    A =  matvitRiesz2rnd(Omega_,n,nu,T); 
elseif strcmp( dist, 'FRiesz' )
    if length(dfs) ~= k+k
        error('length(dfs) wrong.')
    end    
    n = dfs(1 : k); 
    nu = dfs(k + 1 : k + k);
    A =  matvFRieszrnd(Omega_,n,nu,T); 
elseif strcmp( dist, 'FRiesz2' )
    if length(dfs) ~= k+k
        error('length(dfs) wrong.')
    end    
    n = dfs(1 : k); 
    nu = dfs(k + 1 : k + k);
    A =  matvFRiesz2rnd(Omega_,n,nu,T); 
end


end