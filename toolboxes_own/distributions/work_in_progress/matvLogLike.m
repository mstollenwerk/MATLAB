function [ nLogL, logLcontr, varargout] = ...
    matvLogLike(dist, Omega_, dfs, R)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
k = size(Omega_,1);

if strcmp( dist, 'Wish' )
    if length(dfs) ~= 1
        error('length(dfs) wrong.')
    end    
    n = dfs(1);
    like = @(R) matvWishlike(Omega_,n,R);
elseif strcmp( dist, 'iWish' )
    if length(dfs) ~= 1
        error('length(dfs) wrong.')
    end    
    nu = dfs(1);
    like = @(R) matviWishlike(Omega_,nu,R);   
elseif strcmp( dist, 'tWish' )
    if length(dfs) ~= 2
        error('length(dfs) wrong.')
    end    
    n = dfs(1); 
    nu = dfs(2); 
    like = @(R) matvtWishlike(Omega_,n,nu,R);    
elseif strcmp( dist, 'itWish' )
    if length(dfs) ~= 2
        error('length(dfs) wrong.')
    end        
    n = dfs(1); 
    nu = dfs(2);
    like = @(R) matvitWishlike(Omega_,n,nu,R);
elseif strcmp( dist, 'F' )
    if length(dfs) ~= 2
        error('length(dfs) wrong.')
    end    
    n = dfs(1); 
    nu = dfs(2);
    like = @(R) matvFlike(Omega_,n,nu,R);    
elseif strcmp( dist, 'Riesz' )
    if length(dfs) ~= k
        error('length(dfs) wrong.')
    end        
    n = dfs(1 : k);
    like = @(R) matvRieszlike(Omega_,n,R);
elseif strcmp( dist, 'iRiesz2' )
    if length(dfs) ~= k
        error('length(dfs) wrong.')
    end        
    nu = dfs(1 : k); 
    like = @(R) matviRiesz2like(Omega_,nu,R); 
elseif strcmp( dist, 'tRiesz' )
    if length(dfs) ~= k+1
        error('length(dfs) wrong.')
    end        
    n = dfs(1 : k); 
    nu = dfs(k + 1);
    like = @(R) matvtRieszlike(Omega_,n,nu,R);
elseif strcmp( dist, 'itRiesz2' )
    if length(dfs) ~= k+1
        error('length(dfs) wrong.')
    end    
    n = dfs(1); 
    nu = dfs(2 : k + 1);
    like = @(R) matvitRiesz2like(Omega_,n,nu,R); 
elseif strcmp( dist, 'FRiesz' )
    if length(dfs) ~= k+k
        error('length(dfs) wrong.')
    end    
    n = dfs(1 : k); 
    nu = dfs(k + 1 : k + k);
    like = @(R) matvFRieszlike(Omega_,n,nu,R); 
elseif strcmp( dist, 'iFRiesz2' )
    if length(dfs) ~= k+k
        error('length(dfs) wrong.')
    end    
    n = dfs(1 : k); 
    nu = dfs(k + 1 : k + k);
    like = @(R) matviFRiesz2like(Omega_,n,nu,R);
end

if nargout <= 2
    [ nLogL, logLcontr ] = like(R);
elseif nargout == 3    
    [ nLogL, logLcontr, score ] = like(R);
    varargout{1} = score;
elseif nargout == 4    
    [ nLogL, logLcontr, score, hessian ] = like(R);
    varargout{1} = score;
    varargout{2} = hessian;
elseif nargout == 5    
    [ nLogL, logLcontr, score, hessian, param ] = like(R);
    varargout{1} = score;
    varargout{2} = hessian;
    varargout{3} = param;    
elseif nargout == 6
    [ nLogL, logLcontr, score, hessian, param, fisherinfo ] = like(R);
    varargout{1} = score;
    varargout{2} = hessian;
    varargout{3} = param;
    varargout{4} = fisherinfo;
end

end