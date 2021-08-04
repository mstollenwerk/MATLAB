function [ vmat_n, vmat_t ] = matvtWish_maybe_varmat(Sigma, df_w, df_t)
%WISHVARMAT calculates covariance matrix of vector vech(wishmatrix)
%
% USAGE:
%  [VMAT] = WISHVARMAT(SIGMA,DF_W,DF_T)
%
% INPUTS:
%   SIGMA        - K by K array, t-distribution scale matrix.
%   DF_W         - Scalar, Wishart degrees of freedom (dof).
%   DF_T         - Scalar, dof of the underlying t-distribution.
%   
% OUTPUTS:
%   VMAT         - K(K+1)/2 by K(K+1)/2 array, Cov(vech(X)), where X is
%                  t Wishart distributed random variable. 
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
% 19.02.2017

k = size(Sigma,1);
vmat_n=NaN(k*(k+1)/2,k*(k+1)/2);
%% Calcualting normal Wishart CovMat.
cnt1=0;
for i=1:k
    for j=i:k
        cnt1=cnt1+1;
        cnt2=0;
        for l=1:k
            for m=l:k
                cnt2=cnt2+1;
                vmat_n(cnt1,cnt2)=(Sigma(i,l)*Sigma(j,m)+Sigma(i,m)*Sigma(j,l)).*df_w;
            end
        end
    end
end
%% Adjusting for t Wishart CovMat.
vmat_t = df_t^2/(df_t-2)/(df_t-4) * df_w * (vmat_n + 2*df_w/(df_t-2) * vech(Sigma)*vech(Sigma)');
end

