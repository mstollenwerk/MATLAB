clear
clc

Omega_ = rndpd(3);
R(:,:,1) = rndpd(3);
R(:,:,2) = rndpd(3);
n_ = randi([4,12],1);
nu_ = randi([4,12],1);
n = randi([4,12],3,1);
nu = randi([4,12],3,1);
k = 3;

fdiff{1} = finitediff(@(dfs) matvWishlike(Omega_,dfs,R), n_, eps^(1/3) );
fdiff{2} = finitediff(@(dfs) matviWishlike(Omega_,dfs,R), nu_, eps^(1/3) );   
fdiff{3} = finitediff(@(dfs) matvtWishlike(Omega_,dfs(1),dfs(2),R), [n_;nu_], eps^(1/3) );    
fdiff{4} = finitediff(@(dfs) matvitWishlike(Omega_,dfs(1),dfs(2),R), [n_;nu_], eps^(1/3) );
fdiff{5} = finitediff(@(dfs) matvFlike(Omega_,dfs(1),dfs(2),R), [n_;nu_], eps^(1/3) );    
fdiff{6} = finitediff(@(dfs) matvRieszlike(Omega_,dfs,R), n, eps^(1/3) );
fdiff{7} = finitediff(@(dfs) matviRiesz2like(Omega_,dfs,R), nu, eps^(1/3) ); 
fdiff{8} = finitediff(@(dfs) matvtRieszlike(Omega_,dfs(1:k),dfs(k+1),R), [n; nu_], eps^(1/3) );
fdiff{9} = finitediff(@(dfs) matvitRiesz2like(Omega_,dfs(1),dfs(2 : k + 1),R), [n_; nu], eps^(1/3) ); 
fdiff{10} = finitediff(@(dfs) matvFRieszlike(Omega_,dfs(1 : k),dfs(k + 1 : k + k),R), [n; nu], eps^(1/3) ); 
fdiff{11} = finitediff(@(dfs) matviFRiesz2like(Omega_,dfs(1 : k),dfs(k + 1 : k + k),R), [n; nu], eps^(1/3) );

fdiff{12} = finitediff(@(dfs) matvsRieszlike(Omega_,dfs,R), n, eps^(1/3) );

g{1} = matvLogLikeGradientTheta('Wish', repmat(Omega_,1,1,2), repmat(n_,2,1), R);
g{2} = matvLogLikeGradientTheta('iWish', repmat(Omega_,1,1,2), repmat(nu_,2,1), R);
g{3} = matvLogLikeGradientTheta('tWish', repmat(Omega_,1,1,2), [repmat(n_,2,1), repmat(nu_,2,1)], R);
g{4} = matvLogLikeGradientTheta('itWish', repmat(Omega_,1,1,2), [repmat(n_,2,1), repmat(nu_,2,1)], R);
g{5} = matvLogLikeGradientTheta('F', repmat(Omega_,1,1,2), [repmat(n_,2,1), repmat(nu_,2,1)], R);
g{6} = matvLogLikeGradientTheta('Riesz', repmat(Omega_,1,1,2), repmat(n',2,1), R);
g{7} = matvLogLikeGradientTheta('iRiesz2', repmat(Omega_,1,1,2), repmat(nu',2,1), R);
g{8} = matvLogLikeGradientTheta('tRiesz', repmat(Omega_,1,1,2), [repmat(n',2,1), repmat(nu_,2,1)], R);
g{9} = matvLogLikeGradientTheta('itRiesz2', repmat(Omega_,1,1,2), [repmat(n_,2,1), repmat(nu',2,1)], R);
g{10} = matvLogLikeGradientTheta('FRiesz', repmat(Omega_,1,1,2), [repmat(n',2,1), repmat(nu',2,1)], R);
g{11} = matvLogLikeGradientTheta('iFRiesz2', repmat(Omega_,1,1,2), [repmat(n',2,1), repmat(nu',2,1)], R);

g{12} = matvLogLikeGradientTheta('sRiesz', repmat(Omega_,1,1,2), repmat(n',2,1), R);

for ii = 1:12
    err(ii) = sum(abs(fdiff{ii}+sum(g{ii})'));
end
err