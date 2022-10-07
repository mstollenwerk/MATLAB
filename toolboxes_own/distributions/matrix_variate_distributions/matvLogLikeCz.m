function [ nLogL, logLcontr] = matvLogLikeCz(dist, Cz, dfs, X)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

[p,~,T_Cz] = size(Cz);
T_dfs = size(dfs,1);

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
    logLcontr = matvWishlikeCz(Cz,n);
elseif strcmp( dist, 'iWish' )
    if size(dfs,2) ~= 1
        error('dfs must be input in T by size(dfs_t)) format.')
    end    
    nu = dfs;
    logLcontr = matviWishlikeCz(Cz,nu,X);   
elseif strcmp( dist, 'tWish' )
    if size(dfs,2) ~= 2
        error('dfs must be input in T by size(dfs_t)) format.')
    end        
    n = dfs(:,1); 
    nu = dfs(:,2); 
    logLcontr = matvtWishlikeCz(Cz,n,nu);    
elseif strcmp( dist, 'itWish' )
    if size(dfs,2) ~= 2
        error('dfs must be input in T by size(dfs_t)) format.')
    end        
    n = dfs(:,1); 
    nu = dfs(:,2);
    logLcontr = matvitWishlikeCz(Cz,n,nu,X);
elseif strcmp( dist, 'F' )
    if size(dfs,2) ~= 2
        error('dfs must be input in T by size(dfs_t)) format.')
    end        
    n = dfs(:,1); 
    nu = dfs(:,2);
    logLcontr = matvFlikeCz(Cz,n,nu,X);
elseif strcmp( dist, 'Riesz' )
    if size(dfs,2) ~= p
        error('dfs must be input in T by size(dfs_t)) format.')
    end 
    n = dfs;
    logLcontr = matvRieszlikeCz(Cz,n);
elseif strcmp( dist, 'iRiesz2' )
    if size(dfs,2) ~= p
        error('dfs must be input in T by size(dfs_t)) format.')
    end 
    nu = dfs;
    logLcontr = matviRiesz2likeCz(Cz,nu,X); 
elseif strcmp( dist, 'tRiesz' )
    if size(dfs,2) ~= p+1
        error('dfs must be input in T by size(dfs_t)) format.')
    end        
    n = dfs(:,1 : p); 
    nu = dfs(:,p + 1);
    logLcontr = matvtRieszlikeCz(Cz,n,nu);
elseif strcmp( dist, 'itRiesz2' )
    if size(dfs,2) ~= p+1
        error('dfs must be input in T by size(dfs_t)) format.')
    end    
    n = dfs(:,1); 
    nu = dfs(:,2 : p + 1);
    logLcontr = matvitRiesz2likeCz(Cz,n,nu,X); 
elseif strcmp( dist, 'FRiesz' )
    if size(dfs,2) ~= p+p
        error('dfs must be input in T by size(dfs_t)) format.')
    end    
    n = dfs(:,1 : p); 
    nu = dfs(:,p + 1 : p + p);
    logLcontr = matvFRieszlikeCz(Cz,n,nu,X); 
elseif strcmp( dist, 'iFRiesz2' )
    if size(dfs,2) ~= p+p
        error('dfs must be input in T by size(dfs_t)) format.')
    end    
    n = dfs(:,1 : p); 
    nu = dfs(:,p + 1 : p + p);
    logLcontr = matviFRiesz2likeCz(Cz,n,nu,X);
end

nLogL = -sum(logLcontr);

end

function logLcontr = matvWishlikeCz( Cz, n )
p = size(Cz,1);
logLcontr = ...
    - n*p/2*log(2) ...
    - mvgammaln(n/2, p) ...
    + n.*sum(log(diag3d(Cz)),2) ...
    - 1/2*reshape(sum(sum(Cz.^2)),[],1);
end

function logLcontr = matvRieszlikeCz( Cz, n )
logLcontr = ...
    - sum(n,2)/2*log(2) ...
    - lgmvgammaln(n./2) ...
    + loglpwdet3d([],n./2, diag3d(Cz)) ...
    - 1/2*reshape(sum(sum(Cz.^2)),[],1);
end

function logLcontr = matviWishlikeCz( Cz, nu, tr_iZ )
p = size(Cz,1);
logLcontr = ...
    - nu*p/2*log(2) ...
    - mvgammaln(nu/2, p) ...
    - nu.*sum(log(diag3d(Cz)),2) ...
    - 1/2*tr_iZ;
end

function logLcontr = matviRiesz2likeCz( Cz, nu, tr_iZ )
logLcontr = ...
    - sum(nu,2)/2*log(2) ...
    - lgmvgammaln(fliplr(nu)./2) ...
    + loglpwdet3d([],-nu./2, diag3d(Cz)) ...
    - 1/2*tr_iZ;
end

function logLcontr = matvtWishlikeCz( Cz, n, nu )
p = size(Cz,1);
logLcontr = ...
    gammaln( (nu + n*p)/2 ) ...
    - gammaln(nu/2) ...
    - mvgammaln(n/2, p) ...
    - p*n./2.*log(nu) ...
    + n.*sum(log(diag3d(Cz)),2) ...
    - (nu + n*p)/2.*log(1+reshape(sum(sum(Cz.^2)),[],1)./nu);
end

function logLcontr = matvtRieszlikeCz( Cz, n, nu )
logLcontr = ...
    gammaln( (nu + sum(n,2))/2 ) ...
    - gammaln(nu/2) ...
    - lgmvgammaln(n./2) ...
    - sum(n,2)./2.*log(nu) ...
    + loglpwdet3d([],n./2, diag3d(Cz)) ...
    - (nu + sum(n,2))./2.*log(1+reshape(sum(sum(Cz.^2)),[],1)./nu);
end

function logLcontr = matvitWishlikeCz( Cz, n, nu, tr_iZ )
p = size(Cz,1);
logLcontr = ...
    gammaln( (n + nu*p)/2 ) ...
    - gammaln(n/2) ...
    - mvgammaln(nu/2, p) ...
    - p*nu./2.*log(n) ...
    - nu.*sum(log(diag3d(Cz)),2) ...
    - (n + nu*p)/2.*log(1+tr_iZ./n);
end

function logLcontr = matvitRiesz2likeCz( Cz, n, nu, tr_iZ )
logLcontr = ...
    gammaln( (n + sum(nu,2))/2 ) ...
    - gammaln(n/2) ...
    - lgmvgammaln(fliplr(nu)./2) ...
    - sum(nu,2)./2.*log(n) ...
    + loglpwdet3d([],-nu./2, diag3d(Cz)) ...
    - (n + sum(nu,2))./2.*log(1+tr_iZ./n);
end

function logLcontr = matvFlikeCz( Cz, n, nu, logdet_IpZ )
p = size(Cz,1);
logLcontr = ...
    - mvbetaln(n/2, nu/2, p) ...
    + n.*sum(log(diag3d(Cz)),2) ...
    - (n+nu)/2.*logdet_IpZ;
end

function logLcontr = matvFRieszlikeCz( Cz, n, nu, diag3d_chol_IpZ )
logLcontr = ...
    lgmvgammaln(fliplr(n + nu)./2) ...
    - lgmvgammaln(fliplr(nu)./2) ...
    - lgmvgammaln(n./2) ...
    + loglpwdet3d([], n./2, diag3d(Cz)) ...
    + loglpwdet3d([], -(n+nu)./2 , diag3d_chol_IpZ);
end

function logLcontr = matviFRiesz2likeCz( Cz, n, nu, diag3d_cholU_IpZtwisted )
logLcontr = ...
    lgmvgammaln((n + nu)./2) ...
    - lgmvgammaln(fliplr(nu)./2) ...
    - lgmvgammaln(n./2) ...
    + loglpwdet3d([], n./2, diag3d(Cz)) ...
    + logupwdet3d([], -(n+nu)./2 , diag3d_cholU_IpZtwisted);
end

% function logLcontr = matviFRiesz2likeCz( Cz, n, nu, diag3d_chol_InvIpZtwisted )
% logLcontr = ...
%     lgmvgammaln((n + nu)./2) ...
%     - lgmvgammaln(fliplr(nu)./2) ...
%     - lgmvgammaln(n./2) ...
%     + loglpwdet3d([], n./2, diag3d(Cz)) ...
%     + loglpwdet3d([], (n+nu)./2 , diag3d_chol_InvIpZtwisted);
% end

% function logLcontr = matviFRiesz2likeCz( Cz, n, nu, diag3d_chol_InvIpIZ )
% logLcontr = ...
%     lgmvgammaln((n + nu)./2) ...
%     - lgmvgammaln(fliplr(nu)./2) ...
%     - lgmvgammaln(n./2) ...
%     + loglpwdet3d([], -nu./2, diag3d(Cz)) ...
%     + loglpwdet3d([], (n+nu)./2 , diag3d_chol_InvIpIZ);
% end

