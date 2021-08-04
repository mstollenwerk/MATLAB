function [ vmat ] = wishvarmat(Sigma_, df, distr)
%WISHVARMAT calculates covariance matrix of vector vech(wishmatrix)
%
% USAGE:
%  [VMAT] = WISHVARMAT(SIGMA,DF,DISTR)
%
% INPUTS:
%   SIGMA_       - K by K by T array, t-distribution scale matrix.
%   DF           - In case of Distribution=='n': Scalar, Wishart degrees of freedom (dof).
%                  In case of Distribution=='t': Array [df_w df_t], Wishart degrees of freedom (dof) and dof of the underlying t-distribution.                     
%   DISTR        - [Optional] Type of Wishart distribution, one of 'n' (Default) or 't'.
%   
% OUTPUTS:
%   VMAT         - K(K+1)/2 by K(K+1)/2  by T array, Cov(vech(X)), where X is
%                  (t-)Wishart distributed random variable. 
%
% COMMENTS:
%
%  See also 
%
% REFERENCES:
%      [1] Muirhead (1982) - Aspects of Multivariate Statistical Theory.
%      [2] Sutradhar and Ali (1989) - A Generalization of the Wishart
%      Distribution for the Elliptical Model and its Moments for the
%      Multivariate t Model. (Mistake in p.d.f. and CovMat)
%      [3] Joarder (1998) - Some useful Wishart expectations based on the
%      multivariate t model. (Correct p.d.f.)

% Michael Stollenwerk
% michael.stollenwerk@live.com
% 30.08.2019

[N,~,T] = size(Sigma_);
%% Calcualting normal Wishart CovMat.
L = sparse(ELmatrix(N));
K = sparse(Cmat rix(N,N));
eye_ = sparse(eye(N^2));
vmat = NaN(N*(N+1)/2, N*(N+1)/2, T);
ps = parallel.Settings; % Only do parfor if...
ps.Pool.AutoCreate = false; % ...parpool has been created.
parfor tt = 1:T
    vmat(:,:,tt) = 1/df*L*(eye_ + K)*kron(Sigma_(:,:,tt),Sigma_(:,:,tt))*L';
%% Adjusting for t Wishart CovMat.
if strcmpi(distr,'t')
    vmat(:,:,tt) = df(2)^2/(df(2)-2)/(df(2)-4) * df(1) * (vmat(:,:,tt) + 2*df(1)/(df(2)-2) * vech(Sigma_(:,:,tt))*vech(Sigma_(:,:,tt))');
end
end
ps.Pool.AutoCreate = true; % Restore default.
end
