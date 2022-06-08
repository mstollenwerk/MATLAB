function EV = matvEV(dist, Sigma_, dfs)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
[k,~,N] = size(Sigma_);
EV = NaN(k,k,N);

if strcmp( dist, 'Wish' )
    if length(dfs) ~= 1
        error('length(dfs) wrong.')
    end
    n = dfs(1);
    for ii = 1:N
        EV(:,:,ii) = Sigma_(:,:,ii)*n;
    end
elseif strcmp( dist, 'iWish' )
    if length(dfs) ~= 1
        error('length(dfs) wrong.')
    end   
    n = dfs(1);
    for ii = 1:N
        EV(:,:,ii) = Sigma_(:,:,ii)/(n-k-1);
    end    
elseif strcmp( dist, 'tWish' )
    if length(dfs) ~= 2
        error('length(dfs) wrong.')
    end    
    n = dfs(1); 
    nu = dfs(2);
    for ii = 1:N
        EV(:,:,ii) = Sigma_(:,:,ii)*n*nu/(nu-2);
    end
elseif strcmp( dist, 'itWish' )
    if length(dfs) ~= 2
        error('length(dfs) wrong.')
    end    
    n = dfs(1); 
    nu = dfs(2);
    for ii = 1:N
        EV(:,:,ii) = Sigma_(:,:,ii)/(n-k-1);
    end
elseif strcmp( dist, 'F' )
    if length(dfs) ~= 2
        error('length(dfs) wrong.')
    end    
    n = dfs(1); 
    nu = dfs(2);
    for ii = 1:N
        EV(:,:,ii) = Sigma_(:,:,ii)*n*nu/(nu-k-1);
    end
elseif strcmp( dist, 'Riesz' )
    if length(dfs) ~= k
        error('length(dfs) wrong.')
    end    
    n = dfs;
    for ii = 1:N
        L = chol(Sigma_(:,:,ii),'lower');
        EV(:,:,ii) = L*diag(n)*L';
    end
elseif strcmp( dist, 'Riesz2' )
    if length(dfs) ~= k
        error('length(dfs) wrong.')
    end    
    n = dfs;
    for ii = 1:N
        U = cholU(Sigma_(:,:,ii));        
        EV(:,:,ii) = U*diag(n)*U';
    end
elseif strcmp( dist, 'iRiesz' )
    if length(dfs) ~= k
        error('length(dfs) wrong.')
    end    
    n = dfs;  
    for ii = 1:N
        U = cholU(Sigma_(:,:,ii));
        EV(:,:,ii) = U*matviRieszexpmat(n)*U';
    end
elseif strcmp( dist, 'iRiesz2' )
    if length(dfs) ~= k
        error('length(dfs) wrong.')
    end    
    nu = dfs; 
    for ii = 1:N
        L = chol(Sigma_(:,:,ii),'lower');
        EV(:,:,ii) = L*matviRiesz2expmat(nu)*L';
    end
elseif strcmp( dist, 'tRiesz' )
    if length(dfs) ~= k + 1
        error('length(dfs) wrong.')
    end    
    n = dfs(1 : k); 
    nu = dfs(k + 1);
    for ii = 1:N
        L = chol(Sigma_(:,:,ii),'lower');
        EV(:,:,ii) = L*(diag(n)*nu/(nu-2))*L';
    end
elseif strcmp( dist, 'tRiesz2' )
    if length(dfs) ~= k + 1
        error('length(dfs) wrong.')
    end    
    n = dfs(1 : k); 
    nu = dfs(k + 1);
    for ii = 1:N
        U = cholU(Sigma_(:,:,ii)); 
        EV(:,:,ii) = U*(diag(n)*nu/(nu-2))*U';
    end
elseif strcmp( dist, 'itRiesz' )
    if length(dfs) ~= k + 1
        error('length(dfs) wrong.')
    end    
    n = dfs(1 : k); 
    nu = dfs(k + 1);
    for ii = 1:N
        U = cholU(Sigma_(:,:,ii));
        EV(:,:,ii) = U*matviRieszexpmat(n)*U';
    end
elseif strcmp( dist, 'itRiesz2' )
    if length(dfs) ~= k + 1
        error('length(dfs) wrong.')
    end    
    n = dfs(1 : k); 
    nu = dfs(k + 1);
    for ii = 1:N
        L = chol(Sigma_(:,:,ii),'lower');
        EV(:,:,ii) = L*matviRiesz2expmat(n)*L';
    end
elseif strcmp( dist, 'FRiesz' )
    if length(dfs) ~= 2*k
        error('length(dfs) wrong.')
    end    
    n = dfs(1 : k); 
    nu = dfs(k + 1 : k + k); 
    for ii = 1:N
        L = chol(Sigma_(:,:,ii),'lower');
        EV(:,:,ii) = L*matvFRieszexpmat(n,nu)*L';
    end
elseif strcmp( dist, 'FRiesz2' )
    if length(dfs) ~= 2*k
        error('length(dfs) wrong.')
    end    
    n = dfs(1 : k); 
    nu = dfs(k + 1 : k + k); 
    for ii = 1:N
        U = cholU(Sigma_(:,:,ii)); 
        EV(:,:,ii) = U*matvFRies2zexpmat(n,nu)*U';
    end
end

end
