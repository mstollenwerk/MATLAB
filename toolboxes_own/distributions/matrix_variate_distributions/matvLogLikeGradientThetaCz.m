function g = matvLogLikeGradientThetaCz(dist, Cz, dfs, X)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

[p,~,T_Cz] = size(Cz);
T_dfs = size(dfs,1);
if T_Cz < T_dfs
    waring('T_Cz < T_dfs has not been thoroughly tested yet.')
end

if T_dfs ~= 1 && T_Cz ~= 1
    if T_Cz ~= T_dfs
        error('Cz and dfs input sizes are incompatible.')
    end
end

if strcmp( dist, 'Wish' )
    if size(dfs,2) ~= 1
        error('dfs must be input in T by size(dfs_t)) format.')
    end    
    n = dfs;
    g = matvWishlikeGradientThetaCz(Cz,n);
elseif strcmp( dist, 'iWish' )
    if size(dfs,2) ~= 1
        error('dfs must be input in T by size(dfs_t)) format.')
    end    
    nu = dfs;
    g = matviWishlikeGradientThetaCz(Cz,nu);   
elseif strcmp( dist, 'tWish' )
    if size(dfs,2) ~= 2
        error('dfs must be input in T by size(dfs_t)) format.')
    end        
    n = dfs(:,1); 
    nu = dfs(:,2); 
    g = matvtWishlikeGradientThetaCz(Cz,n,nu);    
elseif strcmp( dist, 'itWish' )
    if size(dfs,2) ~= 2
        error('dfs must be input in T by size(dfs_t)) format.')
    end        
    n = dfs(:,1); 
    nu = dfs(:,2);
    g = matvitWishlikeGradientThetaCz(Cz,n,nu,X);
elseif strcmp( dist, 'F' )
    if size(dfs,2) ~= 2
        error('dfs must be input in T by size(dfs_t)) format.')
    end        
    n = dfs(:,1); 
    nu = dfs(:,2);
    g = matvFlikeGradientThetaCz(Cz,n,nu,X);
elseif strcmp( dist, 'Riesz' )
    if size(dfs,2) ~= p
        error('dfs must be input in T by size(dfs_t)) format.')
    end 
    n = dfs;
    g = matvRieszlikeGradientThetaCz(Cz,n);
elseif strcmp( dist, 'iRiesz2' )
    if size(dfs,2) ~= p
        error('dfs must be input in T by size(dfs_t)) format.')
    end 
    nu = dfs;
    g = matviRiesz2likeGradientThetaCz(Cz,nu); 
elseif strcmp( dist, 'tRiesz' )
    if size(dfs,2) ~= p+1
        error('dfs must be input in T by size(dfs_t)) format.')
    end        
    n = dfs(:,1 : p); 
    nu = dfs(:,p + 1);
    g = matvtRieszlikeGradientThetaCz(Cz,n,nu);
elseif strcmp( dist, 'itRiesz2' )
    if size(dfs,2) ~= p+1
        error('dfs must be input in T by size(dfs_t)) format.')
    end    
    n = dfs(:,1); 
    nu = dfs(:,2 : p + 1);
    g = matvitRiesz2likeGradientThetaCz(Cz,n,nu,X); 
elseif strcmp( dist, 'FRiesz' )
    if size(dfs,2) ~= p+p
        error('dfs must be input in T by size(dfs_t)) format.')
    end    
    n = dfs(:,1 : p); 
    nu = dfs(:,p + 1 : p + p);
    g = matvFRieszlikeGradientThetaCz(Cz,n,nu,X); 
elseif strcmp( dist, 'iFRiesz2' )
    if size(dfs,2) ~= p+p
        error('dfs must be input in T by size(dfs_t)) format.')
    end    
    n = dfs(:,1 : p); 
    nu = dfs(:,p + 1 : p + p);
    g = matviFRiesz2likeGradientThetaCz(Cz,n,nu,X);
end

end

function g = matvWishlikeGradientThetaCz( Cz, n )
p = size(Cz,1);
g = - p/2*log(2) ...
    - 1/2*sum(mvpsi_summands_2d(repmat(n./2,1,p)),2) ...
    + sum(log(diag3d(Cz)),2);
end

function g = matvRieszlikeGradientThetaCz( Cz, n )
g = - 1/2*log(2) ...
    - 1/2*mvpsi_summands_2d(n./2) ...
    + log(diag3d(Cz));
end

function g = matviWishlikeGradientThetaCz( Cz, nu )
p = size(Cz,1);
g = - p/2*log(2) ...
    - 1/2*sum(mvpsi_summands_2d(repmat(nu./2,1,p)),2) ... 
    - sum(log(diag3d(Cz)),2);
end

function g = matviRiesz2likeGradientThetaCz( Cz, nu )
g = - 1/2*log(2) ...
    - 1/2*fliplr(mvpsi_summands_2d(fliplr(nu)./2)) ...    
    - log(diag3d(Cz));
end

function g = matvtWishlikeGradientThetaCz( Cz, n, nu )
p = size(Cz,1);
trZ = squeeze(sum(sum(Cz.^2)));

g_n = p/2*psi((nu+p*n)/2) ...
      - 1/2*sum(mvpsi_summands_2d(repmat(n./2,1,p)),2) ...
      - p/2*log(nu) ...
      + sum(log(diag3d(Cz)),2) ...
      - p/2*log(1+trZ./nu);

g_nu = 1/2*psi((nu+p*n)/2) ...
       - 1/2*psi(nu/2) ...
       - 1/2*p*n./nu ...
       - 1/2*log(1+trZ./nu) ...
       + 1/2*(nu+p*n).*(trZ./nu.^2)./(1+trZ./nu);

g = [g_n, g_nu];

end

function g = matvtRieszlikeGradientThetaCz( Cz, n, nu )

trZ = squeeze(sum(sum(Cz.^2)));

g_n = 1/2*psi((nu+sum(n,2))/2) ...
      - 1/2*mvpsi_summands_2d(n./2) ...
      - 1/2*log(nu) ...
      + log(diag3d(Cz)) ...
      - 1/2*log(1+trZ./nu);

g_nu = 1/2*psi((nu+sum(n,2))/2) ...
       - 1/2*psi(nu/2) ...
       - 1/2*sum(n,2)./nu ...
       - 1/2*log(1+trZ./nu) ...
       + 1/2*(nu+sum(n,2))./(1+trZ./nu).*trZ./nu.^2;
   
g = [g_n, g_nu];

end

function g = matvitWishlikeGradientThetaCz( Cz, n, nu, tr_iZ )
p = size(Cz,1);

g_n = 1/2*psi((n+p*nu)/2) ...
      - 1/2*psi(n/2) ...
      - 1/2*p*nu./n ...
      - 1/2*log(1+tr_iZ./n) ...
      + 1/2*(n+p*nu)./(1+tr_iZ./n).*(tr_iZ./n.^2);

g_nu = p/2*psi((n+p*nu)/2) ...
       - 1/2*sum(mvpsi_summands_2d(repmat(nu./2,1,p)),2) ...
       - p/2*log(n) ...
       - sum(log(diag3d(Cz)),2) ...
       - p/2*log(1+tr_iZ./n);

g = [g_n, g_nu];

end

function g = matvitRiesz2likeGradientThetaCz( Cz, n, nu, tr_iZ )
g_n = 1/2*psi((n+sum(nu,2))/2) ...
      - 1/2*psi(n/2) ...
      - 1/2*sum(nu,2)./n ...
      - 1/2*log(1+tr_iZ./n) ...
      + 1/2*(n+sum(nu,2))./(1+tr_iZ./n).*(tr_iZ./n.^2);

g_nu = 1/2*psi((n+sum(nu,2))/2) ...
       - 1/2*fliplr(mvpsi_summands_2d(fliplr(nu)./2)) ...
       - 1/2*log(n) ...
       - log(diag3d(Cz)) ...
       - 1/2*log(1+tr_iZ./n);
   
g = [g_n, g_nu];
end

function g = matvFlikeGradientThetaCz( Cz, n, nu, logdet_IpZ )
p = size(Cz,1);

term1 = 1/2*sum(mvpsi_summands_2d(repmat((n+nu)./2,1,p)),2);
g_n = term1 ...
      - 1/2*sum(mvpsi_summands_2d(repmat(n./2,1,p)),2) ...
      + sum(log(diag3d(Cz)),2) ...
      - 1/2*logdet_IpZ;

g_nu = term1 ...
       - 1/2*sum(mvpsi_summands_2d(repmat(nu./2,1,p)),2) ...
       - 1/2*logdet_IpZ;

g = [g_n, g_nu];
end

function g = matvFRieszlikeGradientThetaCz( Cz, n, nu, diag3d_chol_IpZ )
term1 = 1/2*fliplr(mvpsi_summands_2d(fliplr(n+nu)./2));
g_n = term1 ...
      - 1/2*mvpsi_summands_2d(n./2) ...
      + log(diag3d(Cz)) ...
      - log(diag3d_chol_IpZ);

g_nu = term1 ...
       - 1/2*fliplr(mvpsi_summands_2d(fliplr(nu)./2)) ...
       - log(diag3d_chol_IpZ);

g = [g_n, g_nu];
end

function g = matviFRiesz2likeGradientThetaCz( Cz, n, nu, diag3d_chol_invIpInvZ )
term1 = 1/2*mvpsi_summands_2d((n+nu)./2);
g_n = term1 ...
      - 1/2*mvpsi_summands_2d(n./2) ...
      + log(diag3d_chol_invIpInvZ);

g_nu = term1 ...
   - 1/2*fliplr(mvpsi_summands_2d(fliplr(nu)./2)) ...
   - log(diag3d(Cz)) ...
   + log(diag3d_chol_invIpInvZ);
   
g = [g_n, g_nu];
end