clear
clc

Sigma_ = rndpd(5);
L = chol(Sigma_,'lower');
X(:,:,1) = rndpd(5);
X(:,:,2) = rndpd(5);

n_ = 10;
nu_ = 5;

n = ones(5,1)*n_;
nu = ones(5,1)*nu_;

%%
[ nLogL, logLcontr, Score, ~, param] = matvFRieszlike(Sigma_, n, nu, X);
[ nLogL_, logLcontr_, Score_, ~, param_] = matvFlike(Sigma_/nu_, n_, nu_, X);
nLogL-nLogL_
logLcontr-logLcontr_
Score.Sigma_ - Score_.Sigma_
%%
% [ nLogL, logLcontr, Score, ~, param] = matvFRieszlike(L/matvFRieszexpmat(n,1e25*ones(5,1))*L', n, 1e25*ones(5,1), X);
% [ nLogL_, logLcontr_, Score_, ~, param_] = matvRieszlike(L/diag(n)*L', n, X);
% nLogL-nLogL_
% logLcontr-logLcontr_
% Score.Sigma_ - Score_.Sigma_
%%
[ nLogL, logLcontr, Score, ~, param] = matvRieszlike(Sigma_, n, X);
[ nLogL_, logLcontr_, Score_, ~, param_] = matvWishlike(Sigma_, n_, X);
nLogL-nLogL_
logLcontr-logLcontr_
Score.Sigma_ - Score_.Sigma_
%%
[ nLogL, logLcontr, Score, ~, param] = matviRieszlike(Sigma_, n, X);
[ nLogL_, logLcontr_, Score_, ~, param_] = matviWishlike(Sigma_, n_, X);
nLogL-nLogL_
logLcontr-logLcontr_
Score.Sigma_ - Score_.Sigma_
%%
[ nLogL, logLcontr, Score, ~, param] = matvtRieszlike(Sigma_, n, nu_, X);
[ nLogL_, logLcontr_, Score_, ~, param_] = matvtWishlike(Sigma_, n_, nu_, X);
nLogL-nLogL_
logLcontr-logLcontr_
Score.Sigma_ - Score_.Sigma_
%%
[ nLogL, logLcontr, Score, ~, param] = matvitRieszlike(Sigma_, n, nu_, X);
[ nLogL_, logLcontr_, Score_, ~, param_] = matvitWishlike(Sigma_, n_, nu_, X);
nLogL-nLogL_
logLcontr-logLcontr_
Score.Sigma_ - Score_.Sigma_
%%
[ nLogL, logLcontr, Score, ~, param] = matvtRieszlike(Sigma_, n, nu_, X);
[ nLogL_, logLcontr_, Score_, ~, param_] = matvtWishlike(Sigma_, n_, nu_, X);
nLogL-nLogL_
logLcontr-logLcontr_
Score.Sigma_ - Score_.Sigma_