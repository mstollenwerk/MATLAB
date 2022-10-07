clear
clc

Omega_ = rndpd(3);
R(:,:,1) = rndpd(3);
R(:,:,2) = rndpd(3);
n_ = 8;
nu_ = 9;
n = [ 4 7 5 ]';
nu = [ 7 9 8 ]';
k = 3;

Hfdiff{1} = hessian_2sided(@(dfs) matvWishlike(Omega_,dfs,R), n_ );
Hfdiff{2} = hessian_2sided(@(dfs) matviWishlike(Omega_,dfs,R), nu_ );   
Hfdiff{3} = hessian_2sided(@(dfs) matvtWishlike(Omega_,dfs(1),dfs(2),R), [n_;nu_] );    
Hfdiff{4} = hessian_2sided(@(dfs) matvitWishlike(Omega_,dfs(1),dfs(2),R), [n_;nu_] );
Hfdiff{5} = hessian_2sided(@(dfs) matvFlike(Omega_,dfs(1),dfs(2),R), [n_;nu_] );    
Hfdiff{6} = hessian_2sided(@(dfs) matvRieszlike(Omega_,dfs,R), n );
Hfdiff{7} = hessian_2sided(@(dfs) matviRiesz2like(Omega_,dfs,R), nu ); 
Hfdiff{8} = hessian_2sided(@(dfs) matvtRieszlike(Omega_,dfs(1:k),dfs(k+1),R), [n; nu_] );
Hfdiff{9} = hessian_2sided(@(dfs) matvitRiesz2like(Omega_,dfs(1),dfs(2 : k + 1),R), [n_; nu] ); 
Hfdiff{10} = hessian_2sided(@(dfs) matvFRieszlike(Omega_,dfs(1 : k),dfs(k + 1 : k + k),R), [n; nu] ); 
Hfdiff{11} = hessian_2sided(@(dfs) matviFRiesz2like(Omega_,dfs(1 : k),dfs(k + 1 : k + k),R), [n; nu] );

H{1} = matvLogLikeHessianTheta('Wish', repmat(Omega_,1,1,2), repmat(n_,2,1), R);
H{2} = matvLogLikeHessianTheta('iWish', repmat(Omega_,1,1,2), repmat(nu_,2,1), R);
H{3} = matvLogLikeHessianTheta('tWish', repmat(Omega_,1,1,2), [repmat(n_,2,1), repmat(nu_,2,1)], R);
H{4} = matvLogLikeHessianTheta('itWish', repmat(Omega_,1,1,2), [repmat(n_,2,1), repmat(nu_,2,1)], R);
H{5} = matvLogLikeHessianTheta('F', repmat(Omega_,1,1,2), [repmat(n_,2,1), repmat(nu_,2,1)], R);
H{6} = matvLogLikeHessianTheta('Riesz', repmat(Omega_,1,1,2), repmat(n',2,1), R);
H{7} = matvLogLikeHessianTheta('iRiesz2', repmat(Omega_,1,1,2), repmat(nu',2,1), R);
H{8} = matvLogLikeHessianTheta('tRiesz', repmat(Omega_,1,1,2), [repmat(n',2,1), repmat(nu_,2,1)], R);
H{9} = matvLogLikeHessianTheta('itRiesz2', repmat(Omega_,1,1,2), [repmat(n_,2,1), repmat(nu',2,1)], R);
H{10} = matvLogLikeHessianTheta('FRiesz', repmat(Omega_,1,1,2), [repmat(n',2,1), repmat(nu',2,1)], R);
H{11} = matvLogLikeHessianTheta('iFRiesz2', repmat(Omega_,1,1,2), [repmat(n',2,1), repmat(nu',2,1)], R);

for ii = 1:11
    err(ii) = norm(Hfdiff{ii}+sum(H{ii},3)');
end
err