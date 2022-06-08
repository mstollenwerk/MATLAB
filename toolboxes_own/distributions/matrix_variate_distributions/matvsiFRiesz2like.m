function [ nLogL, logLcontr, score, hessian, param ] = matvsiFRiesz2like( Sigma_, n, nu, X, varargin )
%
% Michael Stollenwerk
% michael.stollenwerk@live.com
% 28.04.2022

A = matviFRiesz2expmat(n,nu);
C = chol(Sigma_, 'lower');
C_Omega = C/sqrtm(A);
Omega_ = C_Omega*C_Omega';

[nLogL, logLcontr, score, hessian, param] = matviFRiesz2like(Omega_, n, nu, X);
param.Sigma_ = Sigma_;

for ii = 1:size(X,3)
    B = inv(inv(Omega_) + inv(X(:,:,ii)));
    C_B = chol(B,'lower');
    
    Nabla = -(C'\trilHalfDiag(C'*tril(C'\diag(n) - C'\A/C*C_B*diag(n+nu)*C_B'/C'))/C);
    score.SigmaNonSym(:,:,ii) = Nabla;

    % Accounting for symmetry of Sigma_:    
    score.Sigma_(:,:,ii) = Nabla+Nabla' - diag(diag(Nabla));
end
        
% p = length(n);
% [~,Gplus] = Dmatrix(p);
% F = ELmatrix(p);
% dvechOmega_dvechSigmat = Gplus*kron(C/A,eye(p))*F'/(Gplus*kron(C,eye(p))*F');
% for ii = 1:size(X,3)
%     score.Sigma_(:,:,ii) = ivech(vech(score.Omega_(:,:,ii))'*dvechOmega_dvechSigmat);
% end

end
