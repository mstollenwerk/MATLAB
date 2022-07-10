function [ nLogL, logLcontr, varargout] = ...
    matvsLogLike(dist, Sigma_, dfs, R)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
k = size(Sigma_,1);

if strcmp( dist, 'Wish' )
    if length(dfs) ~= 1
        error('length(dfs) wrong.')
    end    
    n = dfs(1);
    like = @(R) matvsWishlike(Sigma_,n,R);
elseif strcmp( dist, 'iWish' )
    if length(dfs) ~= 1
        error('length(dfs) wrong.')
    end    
    nu = dfs(1);
    like = @(R) matvsiWishlike(Sigma_,nu,R);   
elseif strcmp( dist, 'tWish' )
    if length(dfs) ~= 2
        error('length(dfs) wrong.')
    end    
    n = dfs(1); 
    nu = dfs(2); 
    like = @(R) matvstWishlike(Sigma_,n,nu,R);    
elseif strcmp( dist, 'itWish' )
    if length(dfs) ~= 2
        error('length(dfs) wrong.')
    end        
    n = dfs(1); 
    nu = dfs(2);
    like = @(R) matvsitWishlike(Sigma_,n,nu,R);
elseif strcmp( dist, 'F' )
    if length(dfs) ~= 2
        error('length(dfs) wrong.')
    end    
    n = dfs(1); 
    nu = dfs(2);
    like = @(R) matvsFlike(Sigma_,n,nu,R);    
elseif strcmp( dist, 'Riesz' )
    if length(dfs) ~= k
        error('length(dfs) wrong.')
    end        
    n = dfs(1 : k);
    like = @(R) matvsRieszlike(Sigma_,n,R);
elseif strcmp( dist, 'iRiesz2' )
    if length(dfs) ~= k
        error('length(dfs) wrong.')
    end        
    nu = dfs(1 : k); 
    like = @(R) matvsiRiesz2like(Sigma_,nu,R); 
elseif strcmp( dist, 'tRiesz' )
    if length(dfs) ~= k+1
        error('length(dfs) wrong.')
    end        
    n = dfs(1 : k); 
    nu = dfs(k + 1);
    like = @(R) matvstRieszlike(Sigma_,n,nu,R);
elseif strcmp( dist, 'itRiesz2' )
    if length(dfs) ~= k+1
        error('length(dfs) wrong.')
    end    
    n = dfs(1); 
    nu = dfs(2 : k + 1);
    like = @(R) matvsitRiesz2like(Sigma_,n,nu,R); 
elseif strcmp( dist, 'FRiesz' )
    if length(dfs) ~= k+k
        error('length(dfs) wrong.')
    end    
    n = dfs(1 : k); 
    nu = dfs(k + 1 : k + k);
    like = @(R) matvsFRieszlike(Sigma_,n,nu,R); 
elseif strcmp( dist, 'iFRiesz2' )
    if length(dfs) ~= k+k
        error('length(dfs) wrong.')
    end    
    n = dfs(1 : k); 
    nu = dfs(k + 1 : k + k);
    like = @(R) matvsiFRiesz2like(Sigma_,n,nu,R);
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