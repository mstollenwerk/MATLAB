function H = matvLogLikeHessianTheta(dist, SigOrOm, dfs, R)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

if size(SigOrOm,3) ~= size(dfs,1) || size(dfs,1) ~= size(R,3)
    error('Parameter inputs have non-matching sizes.')
end
p = size(R,1);

if strcmp( dist, 'Wish' )
    if size(dfs,2) ~= 1
        error('length(dfs) wrong.')
    end    
    n = dfs(:,1);
    H = matvWishlikeHessianTheta(SigOrOm,n,R);
elseif strcmp( dist, 'sWish' )
    if size(dfs,2) ~= 1
        error('length(dfs) wrong.')
    end    
    n = dfs(:,1);
    H = matvsWishlikeHessianTheta(SigOrOm,n,R);    
elseif strcmp( dist, 'iWish' )
    if size(dfs,2) ~= 1
        error('length(dfs) wrong.')
    end    
    nu = dfs(:,1);
    H = matviWishlikeHessianTheta(SigOrOm,nu,R);
elseif strcmp( dist, 'tWish' )
    if size(dfs,2) ~= 2
        error('length(dfs) wrong.')
    end    
    n = dfs(:,1); 
    nu = dfs(:,2); 
    H = matvtWishlikeHessianTheta(SigOrOm,n,nu,R);
elseif strcmp( dist, 'itWish' )
    if size(dfs,2) ~= 2
        error('length(dfs) wrong.')
    end        
    n = dfs(:,1); 
    nu = dfs(:,2);
    H = matvitWishlikeHessianTheta(SigOrOm,n,nu,R);
elseif strcmp( dist, 'F' )
    if size(dfs,2) ~= 2
        error('length(dfs) wrong.')
    end    
    n = dfs(:,1); 
    nu = dfs(:,2);
    H = matvFlikeHessianTheta(SigOrOm,n,nu,R);    
elseif strcmp( dist, 'Riesz' )
    if size(dfs,2) ~= p
        error('length(dfs) wrong.')
    end        
    n = dfs(:,1 : p);
    H = matvRieszlikeHessianTheta(SigOrOm,n,R);    
elseif strcmp( dist, 'iRiesz2' )
    if size(dfs,2) ~= p
        error('length(dfs) wrong.')
    end        
    nu = dfs(:,1 : p); 
    H = matviRiesz2likeHessianTheta(SigOrOm,nu,R);
elseif strcmp( dist, 'tRiesz' )
    if size(dfs,2) ~= p+1
        error('length(dfs) wrong.')
    end        
    n = dfs(:,1 : p); 
    nu = dfs(:,p + 1);
    H = matvtRieszlikeHessianTheta(SigOrOm,n,nu,R);    
elseif strcmp( dist, 'itRiesz2' )
    if size(dfs,2) ~= p+1
        error('length(dfs) wrong.')
    end    
    n = dfs(:,1); 
    nu = dfs(:,2 : p + 1);
    H = matvitRiesz2likeHessianTheta(SigOrOm,n,nu,R);        
elseif strcmp( dist, 'FRiesz' )
    if size(dfs,2) ~= p+p
        error('length(dfs) wrong.')
    end    
    n = dfs(:,1 : p); 
    nu = dfs(:,p + 1 : p + p);
    H = matvFRieszlikeHessianTheta(SigOrOm,n,nu,R);    
elseif strcmp( dist, 'iFRiesz2' )
    if size(dfs,2) ~= p+p
        error('length(dfs) wrong.')
    end    
    n = dfs(:,1 : p); 
    nu = dfs(:,p + 1 : p + p);
    H = matviFRiesz2likeHessianTheta(SigOrOm,n,nu,R);        
else
    error('Please specify valid distribution.')
end

end

%% Theta Hessian functions for individual distributions.
function H = matvWishlikeHessianTheta(Omega_, n, R)
    
    [p,~,T] = size(R);
    H = NaN(1,1,T);    
    for tt = 1:T
        H(tt) = -.25*sum(mvpsi(ones(p,1)*n(tt)/2,1));
    end

end

function H = matvsWishlikeHessianTheta(Sigma_, n, R)
    
    [p,~,T] = size(R);
    H = NaN(1,1,T);    
    for tt = 1:T
        H(tt) = .5*p/n -.25*sum(mvpsi(ones(p,1)*n(tt)/2,1));
    end

end

function H = matviWishlikeHessianTheta(Omega_, nu, R)
    
    [p,~,T] = size(R);
    H = NaN(1,1,T);    
    for tt = 1:T
        H(tt) = -.25*sum(mvpsi(ones(p,1)*nu(tt)/2,1));
    end

end

function H = matvtWishlikeHessianTheta(Omega_, n, nu, R)
    
    [p,~,T] = size(R);
    H = NaN(2,2,T);    
    for tt = 1:T
        trInvOmR = trace(Omega_(:,:,tt)\R(:,:,tt));
        
        d_n_n = .25*( ...
                     p^2*psi(1,(nu(tt)+p*n(tt))/2) ...
                     - sum(mvpsi(ones(p,1)*n(tt)/2,1)) ...
                 );
         
        d_n_nu = .25*p*psi(1,(nu(tt)+p*n(tt))/2) ...
                 - .5*p/(nu(tt)+trInvOmR);
        
        d_nu_n = .5*( ...
                      .5*p*psi(1,(nu(tt)+p*n(tt))/2) ...
                      - p/nu(tt) ...
                      + p*(trInvOmR/nu(tt)^2)/(1+trInvOmR/nu(tt)) ...
                 );
          
        d_nu_nu = .5*( ...
                       .5*psi(1,(nu(tt)+p*n(tt))/2) ...
                       - .5*psi(1,nu(tt)/2) ...
                       + trInvOmR/nu(tt)/(nu(tt)+trInvOmR) ...
                       + (p*n(tt) - trInvOmR)/(nu(tt) + trInvOmR)^2 ...
                 );          
                     
        H(:,:,tt) = [d_n_n, d_n_nu;
                     d_nu_n, d_nu_nu];
    end

end

function H = matvitWishlikeHessianTheta(Omega_, n, nu, R)
    
    [p,~,T] = size(R);
    H = NaN(2,2,T);    
    for tt = 1:T
        trOmInvR = trace(Omega_(:,:,tt)/R(:,:,tt));
        
        d_nu_nu = .25*( ...
                     p^2*psi(1,(n(tt)+p*nu(tt))/2) ...
                     - sum(mvpsi(ones(p,1)*nu(tt)/2,1)) ...
                 );
         
        d_nu_n = .25*p*psi(1,(n(tt)+p*nu(tt))/2) ...
                 - .5*p/(n(tt)+trOmInvR);
        
        d_n_nu = .5*( ...
                      .5*p*psi(1,(n(tt)+p*nu(tt))/2) ...
                      - p/n(tt) ...
                      + p*(trOmInvR/n(tt)^2)/(1+trOmInvR/n(tt)) ...
                 );
          
        d_n_n = .5*( ...
                       .5*psi(1,(n(tt)+p*nu(tt))/2) ...
                       - .5*psi(1,n(tt)/2) ...
                       + trOmInvR/n(tt)/(n(tt)+trOmInvR) ...
                       + (p*nu(tt) - trOmInvR)/(n(tt) + trOmInvR)^2 ...
                 );          
                     
        H(:,:,tt) = [d_n_n, d_n_nu;
                     d_nu_n, d_nu_nu];
                 
    end

end

function H = matvFlikeHessianTheta(Omega_, n, nu, R)

    [p,~,T] = size(R);
    H = NaN(2,2,T);
    for tt = 1:T
        
        term1 = sum(mvpsi(ones(p,1)*(n(tt)+nu(tt))/2,1));
        
        d_n_n = .25*( term1 - sum(mvpsi(ones(p,1)*n(tt)/2,1)) );
         
        d_n_nu = .25*term1;
        
        d_nu_n = .25*term1;
          
        d_nu_nu = .25*( term1 - sum(mvpsi(ones(p,1)*nu(tt)/2,1)) );  
                     
        H(:,:,tt) = [d_n_n, d_n_nu;
                     d_nu_n, d_nu_nu];
        
    end
    
end

function H = matvRieszlikeHessianTheta(Omega_, n, R)

    T = size(R,3);
    H = NaN(size(R,1),size(R,1),T);
    for tt = 1:T

        H(:,:,tt) = -.25*diag( mvpsi(n(tt,:)'/2,1) );
              
    end
    
end

function H = matviRiesz2likeHessianTheta(Omega_, nu, R)

    T = size(R,3);
    H = NaN(size(R,1),size(R,1),T);
    for tt = 1:T

        H(:,:,tt) = -.25*diag( flip(mvpsi(flip(nu(tt,:)')/2,1)) );
              
    end
    
end

function H = matvtRieszlikeHessianTheta(Omega_, n, nu, R)

    [p,~,T] = size(R);
    H = NaN(size(R,1)+1,size(R,1)+1,T);
    for tt = 1:T
        
        trInvOmR = trace(Omega_(:,:,tt)\R(:,:,tt));
        d_n_n = .25*( ...
                          psi(1,(nu(tt)+sum(n(tt,:)))/2) ...
                          - diag(mvpsi(n(tt,:)'/2,1)) ...
                );
         
        d_n_nu = ones(p,1)*( ...
                        .25*psi(1,(nu(tt)+sum(n(tt,:)))/2) ...
                        - .5/(nu(tt) + trInvOmR) ...
                 );         
         
        d_nu_nu = .25*psi(1,(nu(tt)+sum(n(tt,:)))/2) ...
                  - .25*psi(1,nu(tt)/2) ...
                  + .5*trInvOmR/nu(tt)/(nu(tt) + trInvOmR) ...
                  + .5*(sum(n(tt,:))-trInvOmR)/(nu(tt)+trInvOmR)^2;
          
        d_nu_n = ones(1,p)*( ...
                        .25*psi(1,(nu(tt)+sum(n(tt,:)))/2) ...
                        - .5/(nu(tt) + trInvOmR) ...
                 );     
        
        H(:,:,tt) = [d_n_n, d_n_nu;
                     d_nu_n, d_nu_nu];
              
    end
    
end

function H = matvitRiesz2likeHessianTheta(Omega_, n, nu, R)

    [p,~,T] = size(R);
    H = NaN(size(R,1)+1,size(R,1)+1,T);
    for tt = 1:T
        
        trOmInvR = trace(Omega_(:,:,tt)/R(:,:,tt));
        d_n_n = .25*psi(1,(n(tt)+sum(nu(tt,:)))/2) ...
                  - .25*psi(1,n(tt)/2) ...
                  + .5*trOmInvR/n(tt)/(n(tt) + trOmInvR) ...
                  + .5*(sum(nu(tt,:))-trOmInvR)/(n(tt)+trOmInvR)^2;
         
        d_n_nu = ones(1,p)*( ...
                        .25*psi(1,(n(tt)+sum(nu(tt,:)))/2) ...
                        - .5/(n(tt) + trOmInvR) ...
                 );         
         
        d_nu_nu = .25*( ...
                          psi(1,(n(tt)+sum(nu(tt,:)))/2) ...
                          - diag(flip(mvpsi(flip(nu(tt,:))'/2,1))) ...
                );
          
        d_nu_n = ones(p,1)*( ...
                        .25*psi(1,(n(tt)+sum(nu(tt,:)))/2) ...
                        - .5/(n(tt) + trOmInvR) ...
                 );     
        
        H(:,:,tt) = [d_n_n, d_n_nu;
                     d_nu_n, d_nu_nu];
              
    end
    
end

function H = matvFRieszlikeHessianTheta(Omega_, n, nu, R)

    T = size(R,3);
    H = NaN(2*size(R,1),2*size(R,1),T);
    for tt = 1:T
        
        term1 = flip(mvpsi(flip(n(tt,:)'+nu(tt,:)')/2,1));
        
        d_n_n = .25*diag( term1 - mvpsi(n(tt,:)'/2,1) );
        
        d_n_nu = .25*diag( term1 );
        
        d_nu_n = .25*diag( term1 );
        
        d_nu_nu = .25*diag( term1 - flip(mvpsi(flip(nu(tt,:)'/2),1)) );            
          
        H(:,:,tt) = [d_n_n, d_n_nu;
                     d_nu_n, d_nu_nu];
              
    end
    
end

function H = matviFRiesz2likeHessianTheta(Omega_, n, nu, R)

    T = size(R,3);
    H = NaN(2*size(R,1),2*size(R,1),T);
    for tt = 1:T
        
        term1 = mvpsi((n(tt,:)'+nu(tt,:)')/2,1);
        
        d_n_n = .25*diag( term1 - mvpsi(n(tt,:)'/2,1) );
        
        d_n_nu = .25*diag( term1 );
        
        d_nu_n = .25*diag( term1 );
        
        d_nu_nu = .25*diag( term1 - flip(mvpsi(flip(nu(tt,:)'/2),1)) );            
          
        H(:,:,tt) = [d_n_n, d_n_nu;
                     d_nu_n, d_nu_nu];
              
    end
    
end