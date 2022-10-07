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

g{1} =  matvLogLikeGradientTheta('Wish',Omega_,n_,R);
g{2} =  matvLogLikeGradientTheta('iWish',Omega_,nu_,R);
g{3} =  matvLogLikeGradientTheta('tWish',Omega_,[n_,nu_],R);
g{4} =  matvLogLikeGradientTheta('itWish',Omega_,[n_,nu_],R);
g{5} =  matvLogLikeGradientTheta('F',Omega_,[n_,nu_],R);
g{6} =  matvLogLikeGradientTheta('Riesz',Omega_,n,R);
g{7} =  matvLogLikeGradientTheta('iRiesz2',Omega_,nu,R);
g{8} =  matvLogLikeGradientTheta('tRiesz',Omega_,[n,nu_],R);
g{9} =  matvLogLikeGradientTheta('itRiesz2',Omega_,[n_,nu],R);
g{10} =  matvLogLikeGradientTheta('FRiesz',Omega_,[n,nu],R);
g{11} =  matvLogLikeGradientTheta('iFRiesz2',Omega_,[n,nu],R);

g_{1} =  matvLogLikeGradientThetaCz('Wish',Cz,n_,[]);
g_{2} =  matvLogLikeGradientThetaCz('iWish',Cz,nu_,tr_iZ);
g_{3} =  matvLogLikeGradientThetaCz('tWish',Cz,[n_,nu_],[]);
g_{4} =  matvLogLikeGradientThetaCz('itWish',Cz,[n_,nu_],tr_iZ);
g_{5} =  matvLogLikeGradientThetaCz('F',Cz,[n_,nu_],logdet_IpZ);
g_{6} =  matvLogLikeGradientThetaCz('Riesz',Cz,n,[]);
g_{7} =  matvLogLikeGradientThetaCz('iRiesz2',Cz,nu,tr_iZ);
g_{8} =  matvLogLikeGradientThetaCz('tRiesz',Cz,[n,nu_],[]);
g_{9} =  matvLogLikeGradientThetaCz('itRiesz2',Cz,[n_,nu],tr_iZ);
g_{10} =  matvLogLikeGradientThetaCz('FRiesz',Cz,[n,nu],diag3d_chol_IpZ);
g_{11} =  matvLogLikeGradientThetaCz('iFRiesz2',Cz,[n,nu],diag3d_chol_invIpInvZ);

for ii = 1:11
    norm_(ii) = norm(g{ii}-g_{ii});
end
norm_