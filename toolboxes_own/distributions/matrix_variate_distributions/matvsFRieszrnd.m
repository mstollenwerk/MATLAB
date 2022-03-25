function A = matvsFRieszrnd( Sigma_, n, nu, T)
%MATVRIESZRND
%
% USAGE:
%   
%
% INPUTS:
%   
%
% OUTPUTS:
%   
%  See also 
%
% COMMENTS:
% 
% REFERENCES:
%      [1] Stollenwerk (2020)
%
% Michael Stollenwerk
% michael.stollenwerk@live.com
% 24.03.2022
%
% DEPENDENCIES:
%
p = size(Sigma_,1);
p_ = p*(p+1)/2 - p;
C = chol(Sigma_,'lower');

A = NaN(p,p,T);
for tt = 1:T
    
    chi2_l = NaN(p,1);
    for pp = 1:p
        chi2_l(pp) = sqrt(chi2rnd( n(pp) - pp + 1 ));
    end
    Bl = diag(chi2_l);
    Bl(tril(true(p),-1)) = randn(p_,1);
%     A.Riesz(:,:,tt) = Bl*Bl';
    
    chi2_u = NaN(p,1);
    for pp = 1:p
        chi2_u(pp) = sqrt(chi2rnd( nu(pp) + pp - p ));
    end
    Bu = diag(chi2_u);
    Bu(triu(true(p),1)) = randn(p_,1);
%     A.iRiesz2(:,:,tt) = inv(Bu*Bu');
    
    m = diag(matviRiesz2expmat(nu));
    a = diag(matvFRieszexpmat(n,nu));
    
%     Rep1 = C*(diag(sqrt(a))\Bu'\Bl*Bl'/Bu/diag(sqrt(a)))*C';
    Rep1 = diag(sqrt(a))\(Bu'\Bl*Bl'/Bu)/diag(sqrt(a));
    Rep2 = diag(sqrt(m))\(Bu'\(diag(sqrt(n))\(Bl*Bl')/diag(sqrt(n)))/Bu)/diag(sqrt(m));
    
    A(:,:,tt,1) = (Rep1 + Rep1')./2; % ensure symmetry
    A(:,:,tt,2) = (Rep2 + Rep2')./2; % ensure symmetry
% 
%     A.FRiesz(:,:,tt) = Bu'\Bl*Bl'/Bu;
end

end

