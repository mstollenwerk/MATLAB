clear
clc

Omega_(:,:,1) = rndpd(3);
Omega_(:,:,2) = rndpd(3);
R(:,:,1) = rndpd(3);
R(:,:,2) = rndpd(3);
Z(:,:,1) = chol(Omega_(:,:,1),'lower')\R(:,:,1)/chol(Omega_(:,:,1),'lower')';
Z(:,:,2) = chol(Omega_(:,:,2),'lower')\R(:,:,2)/chol(Omega_(:,:,2),'lower')';
Cz(:,:,1) = chol(Z(:,:,1),'lower');
Cz(:,:,2) = chol(Z(:,:,2),'lower');
tr_iZ(1,:) = trace(inv(Z(:,:,1)));
tr_iZ(2,:) = trace(inv(Z(:,:,2)));
logdet_IpZ(1,:) = logdet(eye(3) + Z(:,:,1));
logdet_IpZ(2,:) = logdet(eye(3) + Z(:,:,2));
diag3d_chol_IpZ(1,:) = diag(chol(eye(3) + Z(:,:,1),'lower'));
diag3d_chol_IpZ(2,:) = diag(chol(eye(3) + Z(:,:,2),'lower'));
diag3d_chol_invIpInvZ(1,:) = diag(chol(inv(eye(3) + inv(Z(:,:,1))),'lower'));
diag3d_chol_invIpInvZ(2,:) = diag(chol(inv(eye(3) + inv(Z(:,:,2))),'lower'));
n_ = randi([4,12],2,1);
nu_ = randi([4,12],2,1);
n = randi([4,12],2,3);
nu = randi([4,12],2,3);
k = 3;

H{1} =  matvLogLikeHessianTheta('Wish',Omega_,n_,R);
H{2} =  matvLogLikeHessianTheta('iWish',Omega_,nu_,R);
H{3} =  matvLogLikeHessianTheta('tWish',Omega_,[n_,nu_],R);
H{4} =  matvLogLikeHessianTheta('itWish',Omega_,[n_,nu_],R);
H{5} =  matvLogLikeHessianTheta('F',Omega_,[n_,nu_],R);
H{6} =  matvLogLikeHessianTheta('Riesz',Omega_,n,R);
H{7} =  matvLogLikeHessianTheta('iRiesz2',Omega_,nu,R);
H{8} =  matvLogLikeHessianTheta('tRiesz',Omega_,[n,nu_],R);
H{9} =  matvLogLikeHessianTheta('itRiesz2',Omega_,[n_,nu],R);
H{10} =  matvLogLikeHessianTheta('FRiesz',Omega_,[n,nu],R);
H{11} =  matvLogLikeHessianTheta('iFRiesz2',Omega_,[n,nu],R);

H_{1} =  matvLogLikeHessianThetaCz('Wish',Cz,n_,[]);
H_{2} =  matvLogLikeHessianThetaCz('iWish',Cz,nu_,tr_iZ);
H_{3} =  matvLogLikeHessianThetaCz('tWish',Cz,[n_,nu_],[]);
H_{4} =  matvLogLikeHessianThetaCz('itWish',Cz,[n_,nu_],tr_iZ);
H_{5} =  matvLogLikeHessianThetaCz('F',Cz,[n_,nu_],logdet_IpZ);
H_{6} =  matvLogLikeHessianThetaCz('Riesz',Cz,n,[]);
H_{7} =  matvLogLikeHessianThetaCz('iRiesz2',Cz,nu,tr_iZ);
H_{8} =  matvLogLikeHessianThetaCz('tRiesz',Cz,[n,nu_],[]);
H_{9} =  matvLogLikeHessianThetaCz('itRiesz2',Cz,[n_,nu],tr_iZ);
H_{10} =  matvLogLikeHessianThetaCz('FRiesz',Cz,[n,nu],diag3d_chol_IpZ);
H_{11} =  matvLogLikeHessianThetaCz('iFRiesz2',Cz,[n,nu],diag3d_chol_invIpInvZ);

norm_ = NaN(11,1);
for ii = [1,2,3,4,5,6,7,8,9,10,11]
    norm_(ii) = norm(sum(H{ii},3)-sum(H_{ii},3));
end
norm_