function E = matvExpMat(dist, dfs, p)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

if strcmp( dist, 'Wish' )
    if length(dfs) ~= 1
        error('length(dfs) wrong.')
    end
    n = dfs(1);
    E = diag(repmat(n,p,1));
elseif strcmp( dist, 'iWish' )
    if length(dfs) ~= 1
        error('length(dfs) wrong.')
    end   
    nu = dfs(1);
    E = diag(repmat(1/(nu-p-1),p,1));
elseif strcmp( dist, 'tWish' )
    if length(dfs) ~= 2
        error('length(dfs) wrong.')
    end    
    n = dfs(1); 
    nu = dfs(2);
    E = diag(repmat(n*nu/(nu-2),p,1));
elseif strcmp( dist, 'itWish' )
    if length(dfs) ~= 2
        error('length(dfs) wrong.')
    end    
    n = dfs(1); 
    nu = dfs(2);
    E = diag(repmat(1/(nu-p-1),p,1));
elseif strcmp( dist, 'F' )
    if length(dfs) ~= 2
        error('length(dfs) wrong.')
    end    
    n = dfs(1); 
    nu = dfs(2);
    E = diag(repmat(n/(nu-p-1),p,1));
elseif strcmp( dist, 'Riesz' )
    if length(dfs) ~= p
        error('length(dfs) wrong.')
    end    
    n = dfs;
    E = diag(n);
elseif strcmp( dist, 'Riesz2' )
    if length(dfs) ~= p
        error('length(dfs) wrong.')
    end    
    n = dfs;      
    E = diag(n);
elseif strcmp( dist, 'iRiesz' )
    if length(dfs) ~= p
        error('length(dfs) wrong.')
    end    
    nu = dfs;
    E = matviRieszexpmat(nu);
elseif strcmp( dist, 'iRiesz2' )
    if length(dfs) ~= p
        error('length(dfs) wrong.')
    end    
    nu = dfs; 
    E = matviRiesz2expmat(nu);
elseif strcmp( dist, 'tRiesz' )
    if length(dfs) ~= p + 1
        error('length(dfs) wrong.')
    end    
    n = dfs(1 : p); 
    nu = dfs(p + 1);
    E = diag(n)*nu/(nu-2);
elseif strcmp( dist, 'tRiesz2' )
    if length(dfs) ~= p + 1
        error('length(dfs) wrong.')
    end    
    n = dfs(1 : p); 
    nu = dfs(p + 1);
    E = diag(n)*nu/(nu-2);
elseif strcmp( dist, 'itRiesz' )
    if length(dfs) ~= p + 1
        error('length(dfs) wrong.')
    end    
    n = dfs(1); 
    nu = dfs(2:p+1);
    E = matviRieszexpmat(nu);
elseif strcmp( dist, 'itRiesz2' )
    if length(dfs) ~= p + 1
        error('length(dfs) wrong.')
    end    
    n = dfs(1); 
    nu = dfs(2:p + 1);
    E = matviRiesz2expmat(nu);
elseif strcmp( dist, 'FRiesz' )
    if length(dfs) ~= 2*p
        error('length(dfs) wrong.')
    end    
    n = dfs(1 : p); 
    nu = dfs(p + 1 : p + p); 
    E = matvFRieszexpmat(n,nu);
elseif strcmp( dist, 'FRiesz2' )
    if length(dfs) ~= 2*p
        error('length(dfs) wrong.')
    end    
    n = dfs(1 : p); 
    nu = dfs(p + 1 : p + p);  
    E = matvFRiesz2expmat(n,nu);
elseif strcmp( dist, 'iFRiesz' )
    if length(dfs) ~= 2*p
        error('length(dfs) wrong.')
    end    
    error('not done yet')
elseif strcmp( dist, 'iFRiesz2' )
    if length(dfs) ~= 2*p
        error('length(dfs) wrong.')
    end    
    n = dfs(1 : p); 
    nu = dfs(p + 1 : p + p); 
    E = matviFRiesz2expmat(n,nu);
end

end

