function g = matvFRieszGradient( parString, Omega_, n, nu, X, varargin )
%
% Michael Stollenwerk
% michael.stollenwerk@live.com
% 24.08.2022

[p,~,N] = size(X);
p_ = p*(p+1)/2;
narginchk(5,6);
nargoutchk(0,1);
%% Param
if nargin == 5
    if ~(isempty(Omega_) && isempty(n) && isempty(nu))
        error('Cannot input all_param and any of the parameters individually!')
    end
    [Omega_, ~, n, nu] = matvParamTrans( 'FRiesz', param, p );
end
COm = chol(Omega_,'lower');

g.Omega_ = NaN(N,p_);
g.n = NaN(N,p);
g.nu = NaN(N,p);
g.n_scaled = NaN(N,p);
g.nu_scaled = NaN(N,p);

for ii = 1:N
    
    if any(strcmp("Omega",parString))
        A = X(:,:,ii);

        C_Omplusa = chol(Omega_ + A, 'lower');
        S = COm'\diag(nu)/COm - C_Omplusa'\diag(nu+n)/C_Omplusa;
        S = .5*S;

        % Accounting for symmetry of Omega_:
        S = 2*S - diag(diag(S));

        % Numerical inaccuracies make S not exactly symmetric
        S = (S+S')./2;

        g.Omega_(ii,:) = vech(S);
    end

    if any(strcmp("n",parString))
        g.n(ii,:) = .5*flip(mvpsi(flip(n+nu)/2)) - .5*mvpsi(n/2) + log(diag(chol(A,'lower'))) - log(diag(C_Omplusa));
        g.n_scaled(ii,:) = -g.n(ii,:)'./( .25*flip(mvpsi(flip(n+nu)/2,1)) - .25*mvpsi(n/2,1) );
    end
    if any(strcmp("nu",parString))
        g.nu(ii,:) = .5*flip(mvpsi(flip(n+nu)/2)) - .5*flip(mvpsi(flip(nu/2))) + log(diag(COm)) - log(diag(C_Omplusa));
        g.nu_scaled(ii,:) = -g.nu(ii,:)'./( .25*flip(mvpsi(flip(n+nu)/2,1)) - .25*flip(mvpsi(flip(nu/2),1)) );
    end

end

end
