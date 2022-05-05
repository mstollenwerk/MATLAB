function [ nLogL, logLcontr, score, varargout] = ...
    matvsLogLike(dist, Sigma_, dfs, dta_mat)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
k = size(Sigma_,1);

if strcmp( dist, 'Wish' )
    if length(dfs) ~= 1
        error('length(dfs) wrong.')
    end    
    n = dfs(1);
    [ nLogL, logLcontr, score] = matvsWishlike(Sigma_,n,dta_mat);
elseif strcmp( dist, 'iWish' )
    if length(dfs) ~= 1
        error('length(dfs) wrong.')
    end    
    n = dfs(1);
    [ nLogL, logLcontr, score] = matvsiWishlike(Sigma_,n,dta_mat);   
elseif strcmp( dist, 'tWish' )
    if length(dfs) ~= 2
        error('length(dfs) wrong.')
    end    
    n = dfs(1); 
    nu = dfs(2); 
    [ nLogL, logLcontr, score] = matvstWishlike(Sigma_,n,nu,dta_mat);    
elseif strcmp( dist, 'itWish' )
    if length(dfs) ~= 2
        error('length(dfs) wrong.')
    end        
    n = dfs(1); 
    nu = dfs(2);
    [ nLogL, logLcontr, score] = matvsitWishlike(Sigma_,n,nu,dta_mat);
elseif strcmp( dist, 'F' )
    if length(dfs) ~= 2
        error('length(dfs) wrong.')
    end    
    n = dfs(1); 
    nu = dfs(2);
    [ nLogL, logLcontr, score] = matvsFlike(Sigma_,n,nu,dta_mat);    
elseif strcmp( dist, 'Riesz' )
    if length(dfs) ~= k
        error('length(dfs) wrong.')
    end        
    n = dfs(1 : k);
    [ nLogL, logLcontr, score] = matvsRieszlike(Sigma_,n,dta_mat);    
% elseif strcmp( dist, 'Riesz2' )
%     if length(dfs) ~= k
%         error('length(dfs) wrong.')
%     end        
%     n = dfs(1 : k);
%     [ nLogL, logLcontr, score] = matvsRiesz2like(Sigma_,n,dta_mat);        
% elseif strcmp( dist, 'iRiesz' )
%     if length(dfs) ~= k
%         error('length(dfs) wrong.')
%     end        
%     n = dfs(1 : k);
%     [ nLogL, logLcontr, score] = matvsiRieszlike(Sigma_,n,dta_mat); 
elseif strcmp( dist, 'iRiesz2' )
    if length(dfs) ~= k
        error('length(dfs) wrong.')
    end        
    n = dfs(1 : k); 
    [ nLogL, logLcontr, score] = matvsiRiesz2like(Sigma_,n,dta_mat); 
elseif strcmp( dist, 'tRiesz' )
    if length(dfs) ~= k+1
        error('length(dfs) wrong.')
    end        
    n = dfs(1 : k); 
    nu = dfs(k + 1);
    [ nLogL, logLcontr, score] = matvstRieszlike(Sigma_,n,nu,dta_mat); 
% elseif strcmp( dist, 'tRiesz2' )
%     if length(dfs) ~= k+1
%         error('length(dfs) wrong.')
%     end    
%     n = dfs(1 : k); 
%     nu = dfs(k + 1); 
%     [ nLogL, logLcontr, score] = matvstRiesz2like(Sigma_,n,nu,dta_mat); 
% elseif strcmp( dist, 'itRiesz' )
%     if length(dfs) ~= k+1
%         error('length(dfs) wrong.')
%     end    
%     n = dfs(1 : k); 
%     nu = dfs(k + 1);
%     [ nLogL, logLcontr, score] = matvsitRieszlike(Sigma_,n,nu,dta_mat); 
elseif strcmp( dist, 'itRiesz2' )
    if length(dfs) ~= k+1
        error('length(dfs) wrong.')
    end    
    n = dfs(1 : k); 
    nu = dfs(k + 1);
    [ nLogL, logLcontr, score] = matvsitRiesz2like(Sigma_,n,nu,dta_mat); 
elseif strcmp( dist, 'FRiesz' )
    if length(dfs) ~= k+k
        error('length(dfs) wrong.')
    end    
    n = dfs(1 : k); 
    nu = dfs(k + 1 : k + k);
    [ nLogL, logLcontr, score] = matvsFRieszlike(Sigma_,n,nu,dta_mat); 
elseif strcmp( dist, 'iFRiesz2' )
    if length(dfs) ~= k+k
        error('length(dfs) wrong.')
    end    
    n = dfs(1 : k); 
    nu = dfs(k + 1 : k + k);
    [ nLogL, logLcontr, score] = matvsiFRiesz2like(Sigma_,n,nu,dta_mat);
% elseif strcmp( dist, 'FRiesz2' )
%     if length(dfs) ~= k+k
%         error('length(dfs) wrong.')
%     end    
%     n = dfs(1 : k); 
%     nu = dfs(k + 1 : k + k);
%     [ nLogL, logLcontr, score] = matvsFRiesz2like(Sigma_,n,nu,dta_mat);
end

if nargout >= 4
    param.Sigma_ = Sigma_;
    param.n = n;
    if exist('nu','var')
        param.nu = nu;
    end
    varargout{1} = param;
end

end