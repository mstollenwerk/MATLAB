function [ nLogL, logLcontr, score, hessian, param ] = matvsiFRiesz2like( Sigma_, n, nu, X, varargin )
%
% Michael Stollenwerk
% michael.stollenwerk@live.com
% 28.04.2022

[p,~,N] = size(X);

A = matviFRiesz2expmat(n,nu);
C = chol(Sigma_, 'lower');
C_Omega = C/sqrtm(A);
Omega_ = C_Omega*C_Omega';

[nLogL, logLcontr, score, hessian, param] = matviFRiesz2like(Omega_, n, nu, X);
param.Sigma_ = Sigma_;

for ii = 1:N
    B = inv(inv(Omega_) + inv(X(:,:,ii)));
    C_B = chol(B,'lower');
    
    Nabla = -(C'\trilHalfDiag(C'*tril(C'\diag(n) - C'\A/C*C_B*diag(n+nu)*C_B'/C'))/C);
    score.SigmaNonSym(:,:,ii) = Nabla;

    % Accounting for symmetry of Sigma_:    
    score.Sigma_(:,:,ii) = Nabla+Nabla' - diag(diag(Nabla));
    score.n_originalpdf(ii,:) = .5*(mvpsi((n+nu)/2) - mvpsi(n/2)) ...
                              - log(diag(C_Omega)) + log(diag(C_B));
    score.n_originalpdf_scaled(ii,:) = -score.n_originalpdf(ii,:)'./( .25*( mvpsi((n+nu)/2,1) - mvpsi(n/2,1) ) );
    score.nu_originalpdf(ii,:) = .5*(mvpsi((n+nu)/2) - flip(mvpsi(flip(nu)/2))) ...
                               - log(diag(chol(X(:,:,ii),'lower'))) + log(diag(C_B));
    score.nu_originalpdf_scaled(ii,:) = -score.nu_originalpdf(ii,:)'./( .25*(mvpsi((n+nu)/2,1) - flip(mvpsi(flip(nu)/2,1))) );
                           
    score.n(ii,:) = NaN(1,p);
    score.nu(ii,:) = NaN(1,p);
end
        
% p = length(n);
% [~,Gplus] = Dmatrix(p);
% F = ELmatrix(p);
% dvechOmega_dvechSigmat = Gplus*kron(C/A,eye(p))*F'/(Gplus*kron(C,eye(p))*F');
% for ii = 1:size(X,3)
%     score.Sigma_(:,:,ii) = ivech(vech(score.Omega_(:,:,ii))'*dvechOmega_dvechSigmat);
% end

end
