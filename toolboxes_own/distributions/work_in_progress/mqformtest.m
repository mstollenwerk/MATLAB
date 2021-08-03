clear
clc
rng(1)
%% Setting the Parameters
load('C:\Users\Stollenwerk\OneDrive - bwstaff\Research\minute_data\rc_data\rcov3d_1.mat', 'rc_3d')

Sigma_ = mean(rc_3d(1:2,1:2,:),3);

k = 2;

Psi_ = randn(k);
Psi_ = Psi_*Psi_';

A = randn(k);
A = A*A';

df_ = size(A,1);

%% Random Matrix Generation
mcrep_ = 1e2;
T = 1e3;

[ mqform_rnd ] = mqformrnd( ...
    repmat(Sigma_,1,1,mcrep_*T), ...
    repmat(Psi_,1,1,mcrep_*T), ...
    repmat(A,1,1,mcrep_*T) ...
);
mqform_rnd = reshape(mqform_rnd,k,k,T,mcrep_);

%% Estimation
eSigma_ = NaN(k,k,mcrep_);
ePsi_ = NaN(k,k,mcrep_);
eA = NaN(k,k,mcrep_);

nLogL = NaN(mcrep_,1);

for ii = 1:mcrep_
    [ eSigma_(:,:,ii), ePsi_(:,:,ii), eA(:,:,ii), nLogL(ii) ] = mqformest( mqform_rnd(:,:,:,ii), 100 );
end