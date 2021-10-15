function [ nLogL, varargout] = ...
    matvLogLike(dist, Sigma_, dfs, dta_mat)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
k = size(Sigma_,1);

if nargout == 6
    if strcmp( dist, 'Wish' )
        if length(dfs) ~= 1
            error('length(dfs) wrong.')
        end    
        n = dfs(1);
        [ nLogL, logLcontr, score, hessian, param, fisherinfo] = matvWishlike(Sigma_,n,dta_mat);
    elseif strcmp( dist, 'iWish' )
        if length(dfs) ~= 1
            error('length(dfs) wrong.')
        end    
        n = dfs(1);
        [ nLogL, logLcontr, score, hessian, param, fisherinfo] = matviWishlike(Sigma_,n,dta_mat);   
    elseif strcmp( dist, 'tWish' )
        if length(dfs) ~= 2
            error('length(dfs) wrong.')
        end    
        n = dfs(1); 
        nu = dfs(2); 
        [ nLogL, logLcontr, score, hessian, param, fisherinfo] = matvtWishlike(Sigma_,n,nu,dta_mat);    
    elseif strcmp( dist, 'itWish' )
        if length(dfs) ~= 2
            error('length(dfs) wrong.')
        end        
        n = dfs(1); 
        nu = dfs(2);
        [ nLogL, logLcontr, score, hessian, param, fisherinfo] = matvitWishlike(Sigma_,n,nu,dta_mat);
    elseif strcmp( dist, 'F' )
        if length(dfs) ~= 2
            error('length(dfs) wrong.')
        end    
        n = dfs(1); 
        nu = dfs(2);
        [ nLogL, logLcontr, score, hessian, param, fisherinfo] = matvFlike(Sigma_,n,nu,dta_mat);    
    elseif strcmp( dist, 'Riesz' )
        if length(dfs) ~= k
            error('length(dfs) wrong.')
        end        
        n = dfs(1 : k);
        [ nLogL, logLcontr, score, hessian, param, fisherinfo] = matvRieszlike(Sigma_,n,dta_mat);    
    elseif strcmp( dist, 'Riesz2' )
        if length(dfs) ~= k
            error('length(dfs) wrong.')
        end        
        n = dfs(1 : k);
        [ nLogL, logLcontr, score, hessian, param, fisherinfo] = matvRiesz2like(Sigma_,n,dta_mat);        
    elseif strcmp( dist, 'iRiesz' )
        if length(dfs) ~= k
            error('length(dfs) wrong.')
        end        
        n = dfs(1 : k);
        [ nLogL, logLcontr, score, hessian, param, fisherinfo] = matviRieszlike(Sigma_,n,dta_mat); 
    elseif strcmp( dist, 'iRiesz2' )
        if length(dfs) ~= k
            error('length(dfs) wrong.')
        end        
        n = dfs(1 : k); 
        [ nLogL, logLcontr, score, hessian, param, fisherinfo] = matviRiesz2like(Sigma_,n,dta_mat); 
    elseif strcmp( dist, 'tRiesz' )
        if length(dfs) ~= k+1
            error('length(dfs) wrong.')
        end        
        n = dfs(1 : k); 
        nu = dfs(k + 1);
        [ nLogL, logLcontr, score, hessian, param, fisherinfo] = matvtRieszlike(Sigma_,n,nu,dta_mat); 
    elseif strcmp( dist, 'tRiesz2' )
        if length(dfs) ~= k+1
            error('length(dfs) wrong.')
        end    
        n = dfs(1 : k); 
        nu = dfs(k + 1); 
        [ nLogL, logLcontr, score, hessian, param, fisherinfo] = matvtRiesz2like(Sigma_,n,nu,dta_mat); 
    elseif strcmp( dist, 'itRiesz' )
        if length(dfs) ~= k+1
            error('length(dfs) wrong.')
        end    
        n = dfs(1 : k); 
        nu = dfs(k + 1);
        [ nLogL, logLcontr, score, hessian, param, fisherinfo] = matvitRieszlike(Sigma_,n,nu,dta_mat); 
    elseif strcmp( dist, 'itRiesz2' )
        if length(dfs) ~= k+1
            error('length(dfs) wrong.')
        end    
        n = dfs(1 : k); 
        nu = dfs(k + 1);
        [ nLogL, logLcontr, score, hessian, param, fisherinfo] = matvitRiesz2like(Sigma_,n,nu,dta_mat); 
    elseif strcmp( dist, 'FRiesz' )
        if length(dfs) ~= k+k
            error('length(dfs) wrong.')
        end    
        n = dfs(1 : k); 
        nu = dfs(k + 1 : k + k);
        [ nLogL, logLcontr, score, hessian, param, fisherinfo] = matvFRieszlike(Sigma_,n,nu,dta_mat); 
    elseif strcmp( dist, 'FRiesz2' )
        if length(dfs) ~= k+k
            error('length(dfs) wrong.')
        end    
        n = dfs(1 : k); 
        nu = dfs(k + 1 : k + k);
        [ nLogL, logLcontr, score, hessian, param, fisherinfo] = matvFRiesz2like(Sigma_,n,nu,dta_mat); 
    end
    varargout{1} = logLcontr;
    varargout{2} = score;
    varargout{3} = hessian;
    varargout{5} = param ;
    varargout{6} = fisherinfo;
elseif nargout == 1
    if strcmp( dist, 'Wish' )
        if length(dfs) ~= 1
            error('length(dfs) wrong.')
        end    
        n = dfs(1);
        nLogL = matvWishlike(Sigma_,n,dta_mat);
    elseif strcmp( dist, 'iWish' )
        if length(dfs) ~= 1
            error('length(dfs) wrong.')
        end    
        n = dfs(1);
        nLogL = matviWishlike(Sigma_,n,dta_mat);   
    elseif strcmp( dist, 'tWish' )
        if length(dfs) ~= 2
            error('length(dfs) wrong.')
        end    
        n = dfs(1); 
        nu = dfs(2); 
        nLogL = matvtWishlike(Sigma_,n,nu,dta_mat);    
    elseif strcmp( dist, 'itWish' )
        if length(dfs) ~= 2
            error('length(dfs) wrong.')
        end        
        n = dfs(1); 
        nu = dfs(2);
        nLogL = matvitWishlike(Sigma_,n,nu,dta_mat);
    elseif strcmp( dist, 'F' )
        if length(dfs) ~= 2
            error('length(dfs) wrong.')
        end    
        n = dfs(1); 
        nu = dfs(2);
        nLogL = matvFlike(Sigma_,n,nu,dta_mat);    
    elseif strcmp( dist, 'Riesz' )
        if length(dfs) ~= k
            error('length(dfs) wrong.')
        end        
        n = dfs(1 : k);
        nLogL = matvRieszlike(Sigma_,n,dta_mat);    
    elseif strcmp( dist, 'Riesz2' )
        if length(dfs) ~= k
            error('length(dfs) wrong.')
        end        
        n = dfs(1 : k);
        nLogL = matvRiesz2like(Sigma_,n,dta_mat);        
    elseif strcmp( dist, 'iRiesz' )
        if length(dfs) ~= k
            error('length(dfs) wrong.')
        end        
        n = dfs(1 : k);
        nLogL = matviRieszlike(Sigma_,n,dta_mat); 
    elseif strcmp( dist, 'iRiesz2' )
        if length(dfs) ~= k
            error('length(dfs) wrong.')
        end        
        n = dfs(1 : k); 
        nLogL = matviRiesz2like(Sigma_,n,dta_mat); 
    elseif strcmp( dist, 'tRiesz' )
        if length(dfs) ~= k+1
            error('length(dfs) wrong.')
        end        
        n = dfs(1 : k); 
        nu = dfs(k + 1);
        nLogL = matvtRieszlike(Sigma_,n,nu,dta_mat); 
    elseif strcmp( dist, 'tRiesz2' )
        if length(dfs) ~= k+1
            error('length(dfs) wrong.')
        end    
        n = dfs(1 : k); 
        nu = dfs(k + 1); 
        nLogL = matvtRiesz2like(Sigma_,n,nu,dta_mat); 
    elseif strcmp( dist, 'itRiesz' )
        if length(dfs) ~= k+1
            error('length(dfs) wrong.')
        end    
        n = dfs(1 : k); 
        nu = dfs(k + 1);
        nLogL = matvitRieszlike(Sigma_,n,nu,dta_mat); 
    elseif strcmp( dist, 'itRiesz2' )
        if length(dfs) ~= k+1
            error('length(dfs) wrong.')
        end    
        n = dfs(1 : k); 
        nu = dfs(k + 1);
        nLogL = matvitRiesz2like(Sigma_,n,nu,dta_mat); 
    elseif strcmp( dist, 'FRiesz' )
        if length(dfs) ~= k+k
            error('length(dfs) wrong.')
        end    
        n = dfs(1 : k); 
        nu = dfs(k + 1 : k + k);
        nLogL = matvFRieszlike(Sigma_,n,nu,dta_mat); 
    elseif strcmp( dist, 'FRiesz2' )
        if length(dfs) ~= k+k
            error('length(dfs) wrong.')
        end    
        n = dfs(1 : k); 
        nu = dfs(k + 1 : k + k);
        nLogL = matvFRiesz2like(Sigma_,n,nu,dta_mat); 
    end    
end

end