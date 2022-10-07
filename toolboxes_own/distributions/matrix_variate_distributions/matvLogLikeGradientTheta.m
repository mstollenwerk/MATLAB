function g = matvLogLikeGradientTheta(dist, SigOrOm, dfs, R)
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
    g = matvWishlikeGradientTheta(SigOrOm,n,R);
elseif strcmp( dist, 'sWish' )
    if size(dfs,2) ~= 1
        error('length(dfs) wrong.')
    end    
    n = dfs(:,1);
    g = matvsWishlikeGradientTheta(SigOrOm,n,R);    
elseif strcmp( dist, 'iWish' )
    if size(dfs,2) ~= 1
        error('length(dfs) wrong.')
    end    
    nu = dfs(:,1);
    g = matviWishlikeGradientTheta(SigOrOm,nu,R);
elseif strcmp( dist, 'tWish' )
    if size(dfs,2) ~= 2
        error('length(dfs) wrong.')
    end    
    n = dfs(:,1); 
    nu = dfs(:,2); 
    g = matvtWishlikeGradientTheta(SigOrOm,n,nu,R);
elseif strcmp( dist, 'itWish' )
    if size(dfs,2) ~= 2
        error('length(dfs) wrong.')
    end        
    n = dfs(:,1); 
    nu = dfs(:,2);
    g = matvitWishlikeGradientTheta(SigOrOm,n,nu,R);
elseif strcmp( dist, 'F' )
    if size(dfs,2) ~= 2
        error('length(dfs) wrong.')
    end    
    n = dfs(:,1); 
    nu = dfs(:,2);
    g = matvFlikeGradientTheta(SigOrOm,n,nu,R);    
elseif strcmp( dist, 'Riesz' )
    if size(dfs,2) ~= p
        error('length(dfs) wrong.')
    end        
    n = dfs(:,1 : p);
    g = matvRieszlikeGradientTheta(SigOrOm,n,R); 
elseif strcmp( dist, 'sRiesz' )
    if size(dfs,2) ~= p
        error('length(dfs) wrong.')
    end        
    n = dfs(:,1 : p);
    g = matvsRieszlikeGradientTheta(SigOrOm,n,R);     
elseif strcmp( dist, 'iRiesz2' )
    if size(dfs,2) ~= p
        error('length(dfs) wrong.')
    end        
    nu = dfs(:,1 : p); 
    g = matviRiesz2likeGradientTheta(SigOrOm,nu,R);
elseif strcmp( dist, 'tRiesz' )
    if size(dfs,2) ~= p+1
        error('length(dfs) wrong.')
    end        
    n = dfs(:,1 : p); 
    nu = dfs(:,p + 1);
    g = matvtRieszlikeGradientTheta(SigOrOm,n,nu,R);    
elseif strcmp( dist, 'itRiesz2' )
    if size(dfs,2) ~= p+1
        error('length(dfs) wrong.')
    end    
    n = dfs(:,1); 
    nu = dfs(:,2 : p + 1);
    g = matvitRiesz2likeGradientTheta(SigOrOm,n,nu,R);        
elseif strcmp( dist, 'FRiesz' )
    if size(dfs,2) ~= p+p
        error('length(dfs) wrong.')
    end    
    n = dfs(:,1 : p); 
    nu = dfs(:,p + 1 : p + p);
    g = matvFRieszlikeGradientTheta(SigOrOm,n,nu,R);    
elseif strcmp( dist, 'iFRiesz2' )
    if size(dfs,2) ~= p+p
        error('length(dfs) wrong.')
    end    
    n = dfs(:,1 : p); 
    nu = dfs(:,p + 1 : p + p);
    g = matviFRiesz2likeGradientTheta(SigOrOm,n,nu,R);        
else
    error('Please specify valid distribution.')
end

end

%% Theta gradient functions for individual distributions.
function g = matvWishlikeGradientTheta(Omega_, n, R)
    
    [p,~,T] = size(R);
    g = NaN(T,1);    
    for tt = 1:T
        g(tt) = .5*( ...
                     logdet(R(:,:,tt)) ...
                     - p*log(2)-logdet(Omega_(:,:,tt)) ...
                     - sum(mvpsi(ones(p,1)*n(tt)/2)) ...
               );   
    end

end

function g = matvsWishlikeGradientTheta(Sigma_, n, R)
    
    [p,~,T] = size(R);
    g = NaN(T,1);    
    for tt = 1:T
        g(tt) = .5*( ...
                     p*log(n) ...
                     + p ...
                     + logdet(R(:,:,tt)) ...
                     - p*log(2)-logdet(Sigma_(:,:,tt)) ...
                     - sum(mvpsi(ones(p,1)*n(tt)/2)) ...
                     - trace(Sigma_(:,:,tt)\R(:,:,tt)) ...
               );   
    end

end

function g = matviWishlikeGradientTheta(Omega_, nu, R)
    
    [p,~,T] = size(R);
    g = NaN(T,1);    
    for tt = 1:T
        g(tt) = .5*( ...
                     - p*log(2) ...
                     - sum(mvpsi(ones(p,1)*nu(tt)/2)) ...
                     + logdet(Omega_(:,:,tt)) ...
                     - logdet(R(:,:,tt)) ...
               );  
    end

end

function g = matvtWishlikeGradientTheta(Omega_, n, nu, R)
    
    [p,~,T] = size(R);
    g = NaN(T,2);    
    for tt = 1:T
        trInvOmR = trace(Omega_(:,:,tt)\R(:,:,tt));
        
        g_n = .5*( ...
                   p*psi((nu(tt)+p*n(tt))/2) ...
                   - p*log(nu(tt)) ...
                   - p*log(1+trInvOmR/nu(tt)) ...
                   - sum(mvpsi(ones(p,1)*n(tt)/2)) ...               .                   
                   - logdet(Omega_(:,:,tt)) ...
                   + logdet(R(:,:,tt)) ...
             );
        
        g_nu = .5*( ...
                   psi((nu(tt)+p*n(tt))/2) ...
                   - psi(nu(tt)/2) ...
                   - p*n(tt)/nu(tt) ...
                   - log(1+trInvOmR/nu(tt)) ...
                   + (nu(tt)+p*n(tt))*(trInvOmR/nu(tt)^2)/(1+trInvOmR/nu(tt)) ...
              );
                     
        g(tt,:) = [g_n, g_nu];
    end

end

function g = matvitWishlikeGradientTheta(Omega_, n, nu, R)
    
    [p,~,T] = size(R);
    g = NaN(T,2);    
    for tt = 1:T
        trOmInvR = trace(Omega_(:,:,tt)/R(:,:,tt));
        
        g_n = .5*( ...
                   psi((n(tt)+p*nu(tt))/2) ...
                   - psi(n(tt)/2) ...                   
                   - p*nu(tt)/n(tt) ...
                   - log(1+trOmInvR/n(tt)) ...
                   + (n(tt)+p*nu(tt))*(trOmInvR/n(tt)^2)/(1+trOmInvR/n(tt)) ...
              );
        
        g_nu = .5*( ...
                    p*psi((n(tt)+p*nu(tt))/2) ...
                    - p*log(n(tt)) ...
                    - p*log(1+trOmInvR/n(tt)) ...
                    - sum(mvpsi(ones(p,1)*nu(tt)/2)) ...                    
                    + logdet(Omega_(:,:,tt)) ...
                    - logdet(R(:,:,tt)) ...
              );
          
        g(tt,:) = [g_n, g_nu];
    end

end

function g = matvFlikeGradientTheta(Omega_, n, nu, R)

    [p,~,T] = size(R);
    g = NaN(T,2);
    for tt = 1:T
        
        logdetOmPlusR = logdet(Omega_(:,:,tt) + R(:,:,tt));
        term1 = sum(mvpsi(ones(p,1)*(n(tt)+nu(tt))/2));
        g_n = .5*( ...
                   term1 ...
                   - sum(mvpsi(ones(p,1)*n(tt)/2)) ...
                   + logdet(R(:,:,tt)) ...
                   - logdetOmPlusR ...
             );

        g_nu = .5*( ...
                    term1 ...
                    - sum(mvpsi(ones(p,1)*nu(tt)/2)) ...
                    + logdet(Omega_(:,:,tt)) ...
                    - logdetOmPlusR ...
              );
        
        g(tt,:) = [g_n, g_nu];
        
    end
    
end

function g = matvRieszlikeGradientTheta(Omega_, n, R)

    T = size(R,3);
    g = NaN(T,size(R,1));
    for tt = 1:T

        g(tt,:) = .5*( ...
                       - log(2) ...
                       - mvpsi(n(tt,:)'/2) ...
                  ) ...
                  + log(diag(chol(R(:,:,tt),'lower'))) ...
                  - log(diag(chol(Omega_(:,:,tt),'lower')));
              
    end
    
end

function g = matvsRieszlikeGradientTheta(Sigma_, n, R)

    T = size(R,3);
    g = NaN(T,size(R,1));
    for tt = 1:T
        
        C = chol(Sigma_(:,:,tt),'lower');
        
        g(tt,:) = .5*( ...
                       log(n(tt,:)') ...
                       + 1 ...
                       - log(2) ...
                       - mvpsi(n(tt,:)'/2) ...
                       - diag(C\R(:,:,tt)/C') ...
                  ) ...
                  + log(diag(chol(R(:,:,tt),'lower'))) ...
                  - log(diag(C));
              
    end
    
end

function g = matviRiesz2likeGradientTheta(Omega_, nu, R)

    T = size(R,3);
    g = NaN(T,size(R,1));
    for tt = 1:T

        g(tt,:) = .5*( ...
                       -log(2) ...
                       - flip(mvpsi(flip(nu(tt,:)')/2)) ...
                  ) ...
                  - log(diag(chol(R(:,:,tt),'lower'))) ...
                  + log(diag(chol(Omega_(:,:,tt),'lower')));
              
    end
    
end

function g = matvtRieszlikeGradientTheta(Omega_, n, nu, R)

    T = size(R,3);
    g = NaN(T,size(R,1)+1);
    for tt = 1:T
        
        trInvOmR = trace(Omega_(:,:,tt)\R(:,:,tt));
        g_n = .5*( ...
                   psi((nu(tt)+sum(n(tt,:)))/2) ...
                   - mvpsi(n(tt,:)'/2) ...
                   - log(nu(tt) + trInvOmR) ...
             ) ...
             - log(diag(chol(Omega_(:,:,tt),'lower'))) ...
             + log(diag(chol(R(:,:,tt),'lower')));
         
        g_nu = .5*( ... 
                    psi((nu(tt)+sum(n(tt,:)))/2) ...
                    - psi(nu(tt)/2) ...
                    - log(1 + trInvOmR/nu(tt)) ...
                    - (sum(n(tt,:))-trInvOmR)/(nu(tt)+trInvOmR) ...
              );
        
        g(tt,:) = [g_n', g_nu'];
              
    end
    
end

function g = matvitRiesz2likeGradientTheta(Omega_, n, nu, R)

    T = size(R,3);
    g = NaN(T,size(R,1)+1);
    for tt = 1:T
        
        trOmInvR = trace(Omega_(:,:,tt)/R(:,:,tt));
        g_n = .5*( ...
                   psi((n(tt)+sum(nu(tt,:)))/2) ...
                   - psi(n(tt)/2) ...
                   - sum(nu(tt,:))/n(tt) ...
                   - log(1+trOmInvR/n(tt)) ...
                   + (n(tt)+sum(nu(tt,:)))*(trOmInvR/n(tt)^2)/(1+trOmInvR/n(tt)) ...
             );
         
        g_nu = .5*( ... 
                    psi((n(tt)+sum(nu(tt,:)))/2) ...
                    - flip(mvpsi(flip(nu(tt,:)')/2)) ...
                    - log(n(tt)+trOmInvR) ...
               ) ...
               + log(diag(chol(Omega_(:,:,tt),'lower'))) ...
               - log(diag(chol(R(:,:,tt),'lower')));                      
          
        g(tt,:) = [g_n', g_nu'];
              
    end
    
end

function g = matvFRieszlikeGradientTheta(Omega_, n, nu, R)

    T = size(R,3);
    g = NaN(T,2*size(R,1));
    for tt = 1:T
        
        term1 = .5*flip(mvpsi(flip(n(tt,:)'+nu(tt,:)')/2));
        logdiagcholOmPlusR = log(diag(chol(Omega_(:,:,tt) + R(:,:,tt),'lower')));
        
        g_n = term1 ...
              - .5*mvpsi(n(tt,:)'/2) ...
              + log(diag(chol(R(:,:,tt),'lower'))) ...
              - logdiagcholOmPlusR;
          
        g_nu = term1 ...
               - .5*flip(mvpsi(flip(nu(tt,:)'/2))) ...
               + log(diag(chol(Omega_(:,:,tt),'lower'))) ...
               - logdiagcholOmPlusR;                   
          
        g(tt,:) = [g_n', g_nu'];
              
    end
    
end

function g = matviFRiesz2likeGradientTheta(Omega_, n, nu, R)

    T = size(R,3);
    g = NaN(T,2*size(R,1));
    for tt = 1:T
        
        term1 = .5*mvpsi((n(tt,:)'+nu(tt,:)')/2);
        B = inv(inv(R(:,:,tt)) + inv(Omega_(:,:,tt)));
        logdiagcholB = log(diag(chol(B,'lower')));
        
        g_n = term1 ...
              - .5*mvpsi(n(tt,:)'/2) ...
              - log(diag(chol(Omega_(:,:,tt),'lower'))) ...
              + logdiagcholB;
          
        g_nu = term1 ...
               - .5*flip(mvpsi(flip(nu(tt,:)'/2))) ...
               - log(diag(chol(R(:,:,tt),'lower'))) ...
               + logdiagcholB;                   
          
        g(tt,:) = [g_n', g_nu'];
              
    end
    
end