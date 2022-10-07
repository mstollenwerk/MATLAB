clear
clc

Omega_(:,:,1) = rndpd(3);
Omega_(:,:,2) = rndpd(3);
R(:,:,1) = rndpd(3);
R(:,:,2) = rndpd(3);
Z(:,:,1) = chol(Omega_(:,:,1),'lower')\R(:,:,1)/chol(Omega_(:,:,1),'lower')';
Z(:,:,2) = chol(Omega_(:,:,2),'lower')\R(:,:,2)/chol(Omega_(:,:,2),'lower')';
Z_2(:,:,1) = chol(Omega_(:,:,1),'lower')\R(:,:,1)/chol(Omega_(:,:,1),'lower')';
Z_2(:,:,2) = chol(Omega_(:,:,1),'lower')\R(:,:,2)/chol(Omega_(:,:,1),'lower')';
Cz(:,:,1) = chol(Z(:,:,1),'lower');
Cz(:,:,2) = chol(Z(:,:,2),'lower');
Cz_2(:,:,1) = chol(Z_2(:,:,1),'lower');
Cz_2(:,:,2) = chol(Z_2(:,:,2),'lower');
tr_iZ(1,:) = trace(inv(Z(:,:,1)));
tr_iZ(2,:) = trace(inv(Z(:,:,2)));
logdet_IpZ(1,:) = logdet(eye(3) + Z(:,:,1));
logdet_IpZ(2,:) = logdet(eye(3) + Z(:,:,2));
diag3d_chol_IpZ(1,:) = diag(chol(eye(3) + Z(:,:,1),'lower'));
diag3d_chol_IpZ(2,:) = diag(chol(eye(3) + Z(:,:,2),'lower'));
diag3d_chol_invIpInvZ(1,:) = diag(chol(inv(eye(3) + inv(Z(:,:,1))),'lower'));
diag3d_chol_invIpInvZ(2,:) = diag(chol(inv(eye(3) + inv(Z(:,:,2))),'lower'));
tr_iZ_2(1,:) = trace(inv(Z_2(:,:,1)));
tr_iZ_2(2,:) = trace(inv(Z_2(:,:,2)));
logdet_IpZ_2(1,:) = logdet(eye(3) + Z_2(:,:,1));
logdet_IpZ_2(2,:) = logdet(eye(3) + Z_2(:,:,2));
diag3d_chol_IpZ_2(1,:) = diag(chol(eye(3) + Z_2(:,:,1),'lower'));
diag3d_chol_IpZ_2(2,:) = diag(chol(eye(3) + Z_2(:,:,2),'lower'));
diag3d_chol_invIpInvZ_2(1,:) = diag(chol(inv(eye(3) + inv(Z_2(:,:,1))),'lower'));
diag3d_chol_invIpInvZ_2(2,:) = diag(chol(inv(eye(3) + inv(Z_2(:,:,2))),'lower'));
n_ = randi([4,12],2,1);
nu_ = randi([4,12],2,1);
n = randi([4,12],2,3);
nu = randi([4,12],2,3);
k = 3;

for ii = 1:2
    like_1(1,ii) = matvWishlike(Omega_(:,:,ii),n_(ii),R(:,:,ii));
    like_1(2,ii) = matviWishlike(Omega_(:,:,ii),nu_(ii),R(:,:,ii));
    like_1(3,ii) = matvtWishlike(Omega_(:,:,ii),n_(ii),nu_(ii),R(:,:,ii));
    like_1(4,ii) = matvitWishlike(Omega_(:,:,ii),n_(ii),nu_(ii),R(:,:,ii));
    like_1(5,ii) = matvFlike(Omega_(:,:,ii),n_(ii),nu_(ii),R(:,:,ii));
    like_1(6,ii) = matvRieszlike(Omega_(:,:,ii),n(ii,:)',R(:,:,ii));
    like_1(7,ii) = matviRiesz2like(Omega_(:,:,ii),nu(ii,:)',R(:,:,ii));
    like_1(8,ii) = matvtRieszlike(Omega_(:,:,ii),n(ii,:)',nu_(ii,:)',R(:,:,ii));
    like_1(9,ii) = matvitRiesz2like(Omega_(:,:,ii),n_(ii,:)',nu(ii,:)',R(:,:,ii));
    like_1(10,ii) = matvFRieszlike(Omega_(:,:,ii),n(ii,:)',nu(ii,:)',R(:,:,ii));
    like_1(11,ii) = matviFRiesz2like(Omega_(:,:,ii),n(ii,:)',nu(ii,:)',R(:,:,ii));
end

for ii = 1:2
    like_2(1,ii) = matvWishlike(Omega_(:,:,1),n_(ii),R(:,:,ii));
    like_2(2,ii) = matviWishlike(Omega_(:,:,1),nu_(ii),R(:,:,ii));
    like_2(3,ii) = matvtWishlike(Omega_(:,:,1),n_(ii),nu_(ii),R(:,:,ii));
    like_2(4,ii) = matvitWishlike(Omega_(:,:,1),n_(ii),nu_(ii),R(:,:,ii));
    like_2(5,ii) = matvFlike(Omega_(:,:,1),n_(ii),nu_(ii),R(:,:,ii));
    like_2(6,ii) = matvRieszlike(Omega_(:,:,1),n(ii,:)',R(:,:,ii));
    like_2(7,ii) = matviRiesz2like(Omega_(:,:,1),nu(ii,:)',R(:,:,ii));
    like_2(8,ii) = matvtRieszlike(Omega_(:,:,1),n(ii,:)',nu_(ii,:)',R(:,:,ii));
    like_2(9,ii) = matvitRiesz2like(Omega_(:,:,1),n_(ii,:)',nu(ii,:)',R(:,:,ii));
    like_2(10,ii) = matvFRieszlike(Omega_(:,:,1),n(ii,:)',nu(ii,:)',R(:,:,ii));
    like_2(11,ii) = matviFRiesz2like(Omega_(:,:,1),n(ii,:)',nu(ii,:)',R(:,:,ii));
end

for ii = 1:2
    like_3(1,ii) = matvWishlike(Omega_(:,:,ii),n_(1),R(:,:,ii));
    like_3(2,ii) = matviWishlike(Omega_(:,:,ii),nu_(1),R(:,:,ii));
    like_3(3,ii) = matvtWishlike(Omega_(:,:,ii),n_(1),nu_(1),R(:,:,ii));
    like_3(4,ii) = matvitWishlike(Omega_(:,:,ii),n_(1),nu_(1),R(:,:,ii));
    like_3(5,ii) = matvFlike(Omega_(:,:,ii),n_(1),nu_(1),R(:,:,ii));
    like_3(6,ii) = matvRieszlike(Omega_(:,:,ii),n(1,:)',R(:,:,ii));
    like_3(7,ii) = matviRiesz2like(Omega_(:,:,ii),nu(1,:)',R(:,:,ii));
    like_3(8,ii) = matvtRieszlike(Omega_(:,:,ii),n(1,:)',nu_(1,:)',R(:,:,ii));
    like_3(9,ii) = matvitRiesz2like(Omega_(:,:,ii),n_(1,:)',nu(1,:)',R(:,:,ii));
    like_3(10,ii) = matvFRieszlike(Omega_(:,:,ii),n(1,:)',nu(1,:)',R(:,:,ii));
    like_3(11,ii) = matviFRiesz2like(Omega_(:,:,ii),n(1,:)',nu(1,:)',R(:,:,ii));
end

[~, like_1_(1,:)] =  matvLogLikeCz('Wish',Cz,n_(:),[]);
[~, like_1_(2,:)] =  matvLogLikeCz('iWish',Cz,nu_(:),tr_iZ);
[~, like_1_(3,:)] =  matvLogLikeCz('tWish',Cz,[n_(:),nu_(:)],[]);
[~, like_1_(4,:)] =  matvLogLikeCz('itWish',Cz,[n_(:),nu_(:)],tr_iZ);
[~, like_1_(5,:)] =  matvLogLikeCz('F',Cz,[n_(:),nu_(:)],logdet_IpZ);
[~, like_1_(6,:)] =  matvLogLikeCz('Riesz',Cz,n(:,:),[]);
[~, like_1_(7,:)] =  matvLogLikeCz('iRiesz2',Cz,nu(:,:),tr_iZ);
[~, like_1_(8,:)] =  matvLogLikeCz('tRiesz',Cz,[n(:,:),nu_(:,:)],[]);
[~, like_1_(9,:)] =  matvLogLikeCz('itRiesz2',Cz,[n_(:,:),nu(:,:)],tr_iZ);
[~, like_1_(10,:)] =  matvLogLikeCz('FRiesz',Cz,[n(:,:),nu(:,:)],diag3d_chol_IpZ);
[~, like_1_(11,:)] =  matvLogLikeCz('iFRiesz2',Cz,[n(:,:),nu(:,:)],diag3d_chol_invIpInvZ);

[~, like_2_(1,:)] =  matvLogLikeCz('Wish',Cz_2,n_(:),[]);
[~, like_2_(2,:)] =  matvLogLikeCz('iWish',Cz_2,nu_(:),tr_iZ_2);
[~, like_2_(3,:)] =  matvLogLikeCz('tWish',Cz_2,[n_(:),nu_(:)],[]);
[~, like_2_(4,:)] =  matvLogLikeCz('itWish',Cz_2,[n_(:),nu_(:)],tr_iZ_2);
[~, like_2_(5,:)] =  matvLogLikeCz('F',Cz_2,[n_(:),nu_(:)],logdet_IpZ_2);
[~, like_2_(6,:)] =  matvLogLikeCz('Riesz',Cz_2,n(:,:),[]);
[~, like_2_(7,:)] =  matvLogLikeCz('iRiesz2',Cz_2,nu(:,:),tr_iZ_2);
[~, like_2_(8,:)] =  matvLogLikeCz('tRiesz',Cz_2,[n(:,:),nu_(:,:)],[]);
[~, like_2_(9,:)] =  matvLogLikeCz('itRiesz2',Cz_2,[n_(:,:),nu(:,:)],tr_iZ_2);
[~, like_2_(10,:)] =  matvLogLikeCz('FRiesz',Cz_2,[n(:,:),nu(:,:)],diag3d_chol_IpZ_2);
[~, like_2_(11,:)] =  matvLogLikeCz('iFRiesz2',Cz_2,[n(:,:),nu(:,:)],diag3d_chol_invIpInvZ_2);

[~, like_3_(1,:)] =  matvLogLikeCz('Wish',Cz,n_(1),[]);
[~, like_3_(2,:)] =  matvLogLikeCz('iWish',Cz,nu_(1),tr_iZ);
[~, like_3_(3,:)] =  matvLogLikeCz('tWish',Cz,[n_(1),nu_(1)],[]);
[~, like_3_(4,:)] =  matvLogLikeCz('itWish',Cz,[n_(1),nu_(1)],tr_iZ);
[~, like_3_(5,:)] =  matvLogLikeCz('F',Cz,[n_(1),nu_(1)],logdet_IpZ);
[~, like_3_(6,:)] =  matvLogLikeCz('Riesz',Cz,n(1,:),[]);
[~, like_3_(7,:)] =  matvLogLikeCz('iRiesz2',Cz,nu(1,:),tr_iZ);
[~, like_3_(8,:)] =  matvLogLikeCz('tRiesz',Cz,[n(1,:),nu_(1,:)],[]);
[~, like_3_(9,:)] =  matvLogLikeCz('itRiesz2',Cz,[n_(1,:),nu(1,:)],tr_iZ);
[~, like_3_(10,:)] =  matvLogLikeCz('FRiesz',Cz,[n(1,:),nu(1,:)],diag3d_chol_IpZ);
[~, like_3_(11,:)] =  matvLogLikeCz('iFRiesz2',Cz,[n(1,:),nu(1,:)],diag3d_chol_invIpInvZ);

like_1_ = -like_1_ - repmat(-2*logdet3d(R)',11,1);
like_2_ = -like_2_ - repmat(-2*logdet3d(R)',11,1);
like_3_ = -like_3_ - repmat(-2*logdet3d(R)',11,1);