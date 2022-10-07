function H = matvLogLikeHessianThetaCz(dist, Cz, dfs, X)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

[p,~,T_Cz] = size(Cz);
T_dfs = size(dfs,1);
T = max(T_dfs,T_Cz);

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
    H = matvWishlikeHessianThetaCz(n,p,T);
elseif strcmp( dist, 'iWish' )
    if size(dfs,2) ~= 1
        error('dfs must be input in T by size(dfs_t)) format.')
    end    
    nu = dfs;
    H = matviWishlikeHessianThetaCz(nu,p,T);   
elseif strcmp( dist, 'tWish' )
    if size(dfs,2) ~= 2
        error('dfs must be input in T by size(dfs_t)) format.')
    end        
    n = dfs(:,1); 
    nu = dfs(:,2); 
    H = matvtWishlikeHessianThetaCz(Cz,n,nu);    
elseif strcmp( dist, 'itWish' )
    if size(dfs,2) ~= 2
        error('dfs must be input in T by size(dfs_t)) format.')
    end        
    n = dfs(:,1); 
    nu = dfs(:,2);
    H = matvitWishlikeHessianThetaCz(Cz,n,nu,X);
elseif strcmp( dist, 'F' )
    if size(dfs,2) ~= 2
        error('dfs must be input in T by size(dfs_t)) format.')
    end        
    n = dfs(:,1); 
    nu = dfs(:,2);
    H = matvFlikeHessianThetaCz(n,nu,p);
elseif strcmp( dist, 'Riesz' )
    if size(dfs,2) ~= p
        error('dfs must be input in T by size(dfs_t)) format.')
    end 
    n = dfs;
    H = matvRieszlikeHessianThetaCz(n,T);
elseif strcmp( dist, 'iRiesz2' )
    if size(dfs,2) ~= p
        error('dfs must be input in T by size(dfs_t)) format.')
    end 
    nu = dfs;
    H = matviRiesz2likeHessianThetaCz(nu,T); 
elseif strcmp( dist, 'tRiesz' )
    if size(dfs,2) ~= p+1
        error('dfs must be input in T by size(dfs_t)) format.')
    end        
    n = dfs(:,1 : p); 
    nu = dfs(:,p + 1);
    H = matvtRieszlikeHessianThetaCz(Cz,n,nu);
elseif strcmp( dist, 'itRiesz2' )
    if size(dfs,2) ~= p+1
        error('dfs must be input in T by size(dfs_t)) format.')
    end    
    n = dfs(:,1); 
    nu = dfs(:,2 : p + 1);
    H = matvitRiesz2likeHessianThetaCz(Cz,n,nu,X); 
elseif strcmp( dist, 'FRiesz' )
    if size(dfs,2) ~= p+p
        error('dfs must be input in T by size(dfs_t)) format.')
    end    
    n = dfs(:,1 : p); 
    nu = dfs(:,p + 1 : p + p);
    H = matvFRieszlikeHessianThetaCz(n,nu); 
elseif strcmp( dist, 'iFRiesz2' )
    if size(dfs,2) ~= p+p
        error('dfs must be input in T by size(dfs_t)) format.')
    end    
    n = dfs(:,1 : p); 
    nu = dfs(:,p + 1 : p + p);
    H = matviFRiesz2likeHessianThetaCz(n,nu);
end

end

function H = matvWishlikeHessianThetaCz( n, p, T )
H = - 1/4*sum(mvpsi_summands_2d(repmat(n./2,1,p),1),2);
H = repmat(H,1,1,T);
end

function H = matvRieszlikeHessianThetaCz( n, T )
H = - 1/4*idiag3d(mvpsi_summands_2d(n./2,1));
H = repmat(H,1,1,T);
end

function H = matviWishlikeHessianThetaCz( nu, p, T )
H = - 1/4*sum(mvpsi_summands_2d(repmat(nu./2,1,p),1),2);
H = repmat(H,1,1,T);
end

function H = matviRiesz2likeHessianThetaCz( nu, T )
H = - 1/4*idiag3d(fliplr(mvpsi_summands_2d(fliplr(nu)./2,1)));
H = repmat(H,1,1,T);
end

function H = matvtWishlikeHessianThetaCz( Cz, n, nu )
[p,~,T_Cz] = size(Cz);
T_n = size(n,1);
T = max([T_Cz,T_n]);
H = NaN(2,2,T);

trZ = squeeze(sum(sum(Cz.^2)));

d_n_n = p^2/4*psi(1,(nu+p*n)/2) ...
        - 1/4*sum(mvpsi_summands_2d(repmat(n./2,1,p),1),2);
H(1,1,:) = d_n_n;

d_n_nu = 1/4*p*psi(1,(nu+p*n)/2) ...
         - 1/2*p./(nu+trZ);
H(2,1,:) = d_n_nu;     

d_nu_n = 1/4*p*psi(1,(nu+p*n)/2) ...
         - 1/2*p./nu ...
         + 1/2*p*(trZ./nu.^2)./(1+trZ./nu);
H(1,2,:) = d_nu_n;     

d_nu_nu = 1/4*psi(1,(nu+p*n)/2) ...
          - 1/4*psi(1,nu/2) ...
          + 1/2*trZ./nu./(nu+trZ) ...
          + 1/2*(p*n - trZ)./(nu + trZ).^2; 
H(2,2,:) = d_nu_nu;

end

function H = matvtRieszlikeHessianThetaCz( Cz, n, nu )
[p,~,T_Cz] = size(Cz);
T_n = size(n,1);
T = max([T_Cz,T_n]);
H = NaN(p+1,p+1,T);

trZ = squeeze(sum(sum(Cz.^2)));

d_n_n = 1/4*reshape(repmat(psi(1,(nu+sum(n,2))/2),1,p^2)',p,p,[]) ...
        - 1/4*idiag3d(mvpsi_summands_2d(n./2,1));
if T_n < T_Cz
    d_n_n = repmat(d_n_n,1,1,T);
end
H(1:p,1:p,:) = d_n_n;

d_n_nu = 1/4*psi(1,(nu+sum(n,2))/2) ...
         - 1/2./(nu+trZ);
d_n_nu = repmat(d_n_nu',p,1);
H(1:p,p+1,:) = d_n_nu;     

% d_nu_n = 1/4*p*psi(1,(nu+sum(n,2))/2) ...
%          - 1/2*p./nu ...
%          + 1/2*p*(trZ./nu.^2)./(1+trZ./nu);
% d_nu_n = repmat(d_nu_n,1,p);     
% H(p+1,1:p,:) = d_nu_n';
H(p+1,1:p,:) = d_n_nu; 

d_nu_nu = 1/4*psi(1,(nu+sum(n,2))/2) ...
          - 1/4*psi(1,nu/2) ...
          + 1/2*trZ./nu./(nu+trZ) ...
          + 1/2*(sum(n,2) - trZ)./(nu + trZ).^2; 
H(p+1,p+1,:) = d_nu_nu;
end

function H = matvitWishlikeHessianThetaCz( Cz, n, nu, tr_iZ )
[p,~,T_Cz] = size(Cz);
T_n = size(nu,1);
T = max([T_Cz,T_n]);
H = NaN(2,2,T);

d_n_n = 1/4*psi(1,(n+p*nu)/2) ...
        - 1/4*psi(1,n/2) ...
        + 1/2*tr_iZ./n./(n+tr_iZ) ...
        + 1/2*(p*nu - tr_iZ)./(n + tr_iZ).^2; 
H(1,1,:) = d_n_n;

d_n_nu = 1/4*p*psi(1,(n+p*nu)/2) ...
         - 1/2*p./n ...
         + 1/2*p*(tr_iZ./n.^2)./(1+tr_iZ./n);
H(2,1,:) = d_n_nu;  

% d_nu_n = 1/4*p*psi(1,(n+p*nu)/2) ...
%          - 1/2*p./(n+tr_iZ);
% H(1,2,:) = d_nu_n;  
H(1,2,:) = d_n_nu;

d_nu_nu = p^2/4*psi(1,(n+p*nu)/2) ...
          - 1/4*sum(mvpsi_summands_2d(repmat(nu./2,1,p),1),2);
H(2,2,:) = d_nu_nu;

end

function H = matvitRiesz2likeHessianThetaCz( Cz, n, nu, tr_iZ )
[p,~,T_Cz] = size(Cz);
T_n = size(n,1);
T = max([T_Cz,T_n]);
H = NaN(p+1,p+1,T);

d_n_n = 1/4*psi(1,(n+sum(nu,2))/2) ...
          - 1/4*psi(1,n/2) ...
          + 1/2*tr_iZ./n./(n+tr_iZ) ...
          + 1/2*(sum(nu,2) - tr_iZ)./(n + tr_iZ).^2; 
H(1,1,:) = d_n_n;

d_nu_n = 1/4*psi(1,(n+sum(nu,2))/2) ...
         - 1/2./(n+tr_iZ);
d_nu_n = repmat(d_nu_n',p,1);
H(1,2:p+1,:) = d_nu_n;     

% d_n_nu = 1/4*p*psi(1,(nu+sum(n,2))/2) ...
%          - 1/2*p./nu ...
%          + 1/2*p*(trZ./nu.^2)./(1+trZ./nu);
% d_n_nu = repmat(d_n_nu,1,p);     
% H(1:p,p+1,:) = d_n_nu';
H(2:p+1,1,:) = d_nu_n; 

d_nu_nu = 1/4*reshape(repmat(psi(1,(n+sum(nu,2))/2),1,p^2)',p,p,[]) ...
        - 1/4*idiag3d(fliplr(mvpsi_summands_2d(fliplr(nu)./2,1)));
if T_n < T_Cz
    d_nu_nu = repmat(d_nu_nu,1,1,T);
end
H(2:p+1,2:p+1,:) = d_nu_nu;

end

function H = matvFlikeHessianThetaCz( n, nu, p )
term1 = 1/4*sum(mvpsi_summands_2d(repmat((n+nu)./2,1,p),1),2);
d_n_n = term1 ...
         - 1/4*sum(mvpsi_summands_2d(repmat(n./2,1,p),1),2);

d_n_nu = term1;

d_nu_n = term1;

d_nu_nu = term1 ...
          - 1/4*sum(mvpsi_summands_2d(repmat(nu./2,1,p),1),2);

H = reshape([d_n_n, d_n_nu, d_nu_n, d_nu_nu]',2,2,[]);
if size(H,3) < size(n,1)
    H = repmat(H,1,1,size(n,1));
end
end

function H = matvFRieszlikeHessianThetaCz( n, nu )
term1 = 1/4*fliplr(mvpsi_summands_2d(fliplr(n+nu)./2,1));
d_n_n = idiag3d(term1) ...
         - 1/4*idiag3d(mvpsi_summands_2d(n./2,1));

d_n_nu = idiag3d(term1);

d_nu_n = idiag3d(term1);

d_nu_nu = idiag3d(term1) ...
          - 1/4*idiag3d(fliplr(mvpsi_summands_2d(fliplr(nu)./2,1)));
      
BlockRow1 = cat(2,d_n_n,d_n_nu);
BlockRow2 = cat(2,d_nu_n,d_nu_nu);
H = cat(1,BlockRow1,BlockRow2);

if size(H,3) < size(n,1)
    H = repmat(H,1,1,size(n,1));
end
end

function H = matviFRiesz2likeHessianThetaCz( n, nu )
term1 = 1/4*mvpsi_summands_2d((n+nu)./2,1);
d_n_n = idiag3d(term1) ...
         - 1/4*idiag3d(mvpsi_summands_2d(n./2,1));

d_n_nu = idiag3d(term1);

d_nu_n = idiag3d(term1);

d_nu_nu = idiag3d(term1) ...
          - 1/4*idiag3d(fliplr(mvpsi_summands_2d(fliplr(nu)./2,1)));
      
BlockRow1 = cat(2,d_n_n,d_n_nu);
BlockRow2 = cat(2,d_nu_n,d_nu_nu);
H = cat(1,BlockRow1,BlockRow2);

if size(H,3) < size(n,1)
    H = repmat(H,1,1,size(n,1));
end
end