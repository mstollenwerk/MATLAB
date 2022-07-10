function Sigma_ = matvEV(dist, Omega_, dfs)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
[k,~,N] = size(Omega_);
Sigma_ = NaN(k,k,N);

if strcmp( dist, 'Wish' )
    if length(dfs) ~= 1
        error('length(dfs) wrong.')
    end
    n = dfs(1);
    for ii = 1:N
        Sigma_(:,:,ii) = Omega_(:,:,ii)*n;
    end
elseif strcmp( dist, 'iWish' )
    if length(dfs) ~= 1
        error('length(dfs) wrong.')
    end   
    nu = dfs(1);
    for ii = 1:N
        Sigma_(:,:,ii) = Omega_(:,:,ii)/(nu-k-1);
    end    
elseif strcmp( dist, 'tWish' )
    if length(dfs) ~= 2
        error('length(dfs) wrong.')
    end    
    n = dfs(1); 
    nu = dfs(2);
    for ii = 1:N
        Sigma_(:,:,ii) = Omega_(:,:,ii)*n*nu/(nu-2);
    end
elseif strcmp( dist, 'itWish' )
    if length(dfs) ~= 2
        error('length(dfs) wrong.')
    end    
    n = dfs(1); 
    nu = dfs(2);
    for ii = 1:N
        Sigma_(:,:,ii) = Omega_(:,:,ii)/(nu-k-1);
    end
elseif strcmp( dist, 'F' )
    if length(dfs) ~= 2
        error('length(dfs) wrong.')
    end    
    n = dfs(1); 
    nu = dfs(2);
    for ii = 1:N
        Sigma_(:,:,ii) = Omega_(:,:,ii)*n/(nu-k-1);
    end
elseif strcmp( dist, 'Riesz' )
    if length(dfs) ~= k
        error('length(dfs) wrong.')
    end    
    n = dfs;
    for ii = 1:N
        COm = chol(Omega_(:,:,ii),'lower');
        Sigma_(:,:,ii) = COm*diag(n)*COm';
    end
elseif strcmp( dist, 'Riesz2' )
    if length(dfs) ~= k
        error('length(dfs) wrong.')
    end    
    n = dfs;
    for ii = 1:N
        U = cholU(Omega_(:,:,ii));        
        Sigma_(:,:,ii) = U*diag(n)*U';
    end
elseif strcmp( dist, 'iRiesz' )
    if length(dfs) ~= k
        error('length(dfs) wrong.')
    end    
    nu = dfs;  
    for ii = 1:N
        U = cholU(Omega_(:,:,ii));
        Sigma_(:,:,ii) = U*matviRieszexpmat(nu)*U';
    end
elseif strcmp( dist, 'iRiesz2' )
    if length(dfs) ~= k
        error('length(dfs) wrong.')
    end    
    nu = dfs; 
    for ii = 1:N
        COm = chol(Omega_(:,:,ii),'lower');
        Sigma_(:,:,ii) = COm*matviRiesz2expmat(nu)*COm';
    end
elseif strcmp( dist, 'tRiesz' )
    if length(dfs) ~= k + 1
        error('length(dfs) wrong.')
    end    
    n = dfs(1 : k); 
    nu = dfs(k + 1);
    for ii = 1:N
        COm = chol(Omega_(:,:,ii),'lower');
        Sigma_(:,:,ii) = COm*(diag(n)*nu/(nu-2))*COm';
    end
elseif strcmp( dist, 'tRiesz2' )
    if length(dfs) ~= k + 1
        error('length(dfs) wrong.')
    end    
    n = dfs(1 : k); 
    nu = dfs(k + 1);
    for ii = 1:N
        U = cholU(Omega_(:,:,ii)); 
        Sigma_(:,:,ii) = U*(diag(n)*nu/(nu-2))*U';
    end
elseif strcmp( dist, 'itRiesz' )
    if length(dfs) ~= k + 1
        error('length(dfs) wrong.')
    end    
    n = dfs(1); 
    nu = dfs(2:k+1);
    for ii = 1:N
        U = cholU(Omega_(:,:,ii));
        Sigma_(:,:,ii) = U*matviRieszexpmat(nu)*U';
    end
elseif strcmp( dist, 'itRiesz2' )
    if length(dfs) ~= k + 1
        error('length(dfs) wrong.')
    end    
    n = dfs(1); 
    nu = dfs(2:k+1);
    for ii = 1:N
        COm = chol(Omega_(:,:,ii),'lower');
        Sigma_(:,:,ii) = COm*matviRiesz2expmat(nu)*COm';
    end
elseif strcmp( dist, 'FRiesz' )
    if length(dfs) ~= 2*k
        error('length(dfs) wrong.')
    end    
    n = dfs(1 : k); 
    nu = dfs(k + 1 : k + k); 
    for ii = 1:N
        COm = chol(Omega_(:,:,ii),'lower');
        Sigma_(:,:,ii) = COm*matvFRieszexpmat(n,nu)*COm';
    end
elseif strcmp( dist, 'FRiesz2' )
    if length(dfs) ~= 2*k
        error('length(dfs) wrong.')
    end    
    n = dfs(1 : k); 
    nu = dfs(k + 1 : k + k); 
    for ii = 1:N
        U = cholU(Omega_(:,:,ii)); 
        Sigma_(:,:,ii) = U*matvFRiesz2expmat(n,nu)*U';
    end
elseif strcmp( dist, 'iFRiesz' )
    if length(dfs) ~= 2*k
        error('length(dfs) wrong.')
    end    
    error('not done yet')
elseif strcmp( dist, 'iFRiesz2' )
    if length(dfs) ~= 2*k
        error('length(dfs) wrong.')
    end    
    n = dfs(1 : k); 
    nu = dfs(k + 1 : k + k); 
    for ii = 1:N
        COm = chol(Omega_(:,:,ii),'lower'); 
        Sigma_(:,:,ii) = COm*matviFRiesz2expmat(n,nu)*COm';
    end    
end

end

