function Omega_ = matvStandardize(dist, Sigma_or_R, dfs)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
[k,~,N] = size(Sigma_or_R);
Omega_ = NaN(k,k,N);

if strcmp( dist, 'Wish' )
    if length(dfs) ~= 1
        error('length(dfs) wrong.')
    end
    n = dfs(1);
    for ii = 1:N
        Omega_(:,:,ii) = Sigma_or_R(:,:,ii)/n;
    end
elseif strcmp( dist, 'iWish' )
    if length(dfs) ~= 1
        error('length(dfs) wrong.')
    end   
    nu = dfs(1);
    for ii = 1:N
        Omega_(:,:,ii) = Sigma_or_R(:,:,ii)*(nu-k-1);
    end    
elseif strcmp( dist, 'tWish' )
    if length(dfs) ~= 2
        error('length(dfs) wrong.')
    end    
    n = dfs(1); 
    nu = dfs(2);
    for ii = 1:N
        Omega_(:,:,ii) = Sigma_or_R(:,:,ii)/n/nu*(nu-2);
    end
elseif strcmp( dist, 'itWish' )
    if length(dfs) ~= 2
        error('length(dfs) wrong.')
    end    
    n = dfs(1); 
    nu = dfs(2);
    for ii = 1:N
        Omega_(:,:,ii) = Sigma_or_R(:,:,ii)*(nu-k-1);
    end
elseif strcmp( dist, 'F' )
    if length(dfs) ~= 2
        error('length(dfs) wrong.')
    end    
    n = dfs(1); 
    nu = dfs(2);
    for ii = 1:N
        Omega_(:,:,ii) = Sigma_or_R(:,:,ii)/n*(nu-k-1);
    end
elseif strcmp( dist, 'Riesz' )
    if length(dfs) ~= k
        error('length(dfs) wrong.')
    end    
    n = dfs;
    for ii = 1:N
        C = chol(Sigma_or_R(:,:,ii),'lower');
        Omega_(:,:,ii) = C/diag(n)*C';
    end
elseif strcmp( dist, 'Riesz2' )
    if length(dfs) ~= k
        error('length(dfs) wrong.')
    end    
    n = dfs;
    for ii = 1:N
        U = cholU(Sigma_or_R(:,:,ii));        
        Omega_(:,:,ii) = U/diag(n)*U';
    end
elseif strcmp( dist, 'iRiesz' )
    if length(dfs) ~= k
        error('length(dfs) wrong.')
    end    
    nu = dfs;  
    for ii = 1:N
        U = cholU(Sigma_or_R(:,:,ii));
        Omega_(:,:,ii) = U/matviRieszexpmat(nu)*U';
    end
elseif strcmp( dist, 'iRiesz2' )
    if length(dfs) ~= k
        error('length(dfs) wrong.')
    end    
    nu = dfs; 
    for ii = 1:N
        C = chol(Sigma_or_R(:,:,ii),'lower');
        Omega_(:,:,ii) = C/matviRiesz2expmat(nu)*C';
    end
elseif strcmp( dist, 'tRiesz' )
    if length(dfs) ~= k + 1
        error('length(dfs) wrong.')
    end    
    n = dfs(1 : k); 
    nu = dfs(k + 1);
    for ii = 1:N
        C = chol(Sigma_or_R(:,:,ii),'lower');
        Omega_(:,:,ii) = C/(diag(n)*nu/(nu-2))*C';
    end
elseif strcmp( dist, 'tRiesz2' )
    if length(dfs) ~= k + 1
        error('length(dfs) wrong.')
    end    
    n = dfs(1 : k); 
    nu = dfs(k + 1);
    for ii = 1:N
        U = cholU(Sigma_or_R(:,:,ii)); 
        Omega_(:,:,ii) = U/(diag(n)*nu/(nu-2))*U';
    end
elseif strcmp( dist, 'itRiesz' )
    if length(dfs) ~= k + 1
        error('length(dfs) wrong.')
    end    
    n = dfs(1); 
    nu = dfs(2:k+1);
    for ii = 1:N
        U = cholU(Sigma_or_R(:,:,ii));
        Omega_(:,:,ii) = U/matviRieszexpmat(nu)*U';
    end
elseif strcmp( dist, 'itRiesz2' )
    if length(dfs) ~= k + 1
        error('length(dfs) wrong.')
    end    
    n = dfs(1); 
    nu = dfs(2:k + 1);
    for ii = 1:N
        C = chol(Sigma_or_R(:,:,ii),'lower');
        Omega_(:,:,ii) = C/matviRiesz2expmat(nu)*C';
    end
elseif strcmp( dist, 'FRiesz' )
    if length(dfs) ~= 2*k
        error('length(dfs) wrong.')
    end    
    n = dfs(1 : k); 
    nu = dfs(k + 1 : k + k); 
    for ii = 1:N
        C = chol(Sigma_or_R(:,:,ii),'lower');
        Omega_(:,:,ii) = C/matvFRieszexpmat(n,nu)*C';
    end
elseif strcmp( dist, 'FRiesz2' )
    if length(dfs) ~= 2*k
        error('length(dfs) wrong.')
    end    
    n = dfs(1 : k); 
    nu = dfs(k + 1 : k + k); 
    for ii = 1:N
        U = cholU(Sigma_or_R(:,:,ii)); 
        Omega_(:,:,ii) = U/matvFRiesz2expmat(n,nu)*U';
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
        C = chol(Sigma_or_R(:,:,ii),'lower'); 
        Omega_(:,:,ii) = C/matviFRiesz2expmat(n,nu)*C';
    end
end

end

