function [ X, Sigma_, ScaledScore ] = gas_scalar_BEKK_rnd( intrcpt, scoreparam, garchparam, n, nu, T, dist )
%
% Michael Stollenwerk
% michael.stollenwerk@live.com
% 11.08.2020

k = size(intrcpt,1);
p = length(scoreparam);
q = length(garchparam);
%% Data Storage
X = NaN(k,k,T);
Sigma_ = NaN(k,k,T);
ScaledScore = NaN(k,k,T);
%% Recursion
% Initialize recursion at unconditional mean (stationarity assumed).
ini_Sigma = intrcpt./(1 - sum(garchparam));
for tt=1:T
    % Parameter Recursion
    Sigma_(:,:,tt) = intrcpt;
    for jj = 1:p
        if (tt-jj) <= 0
            Sigma_(:,:,tt) = Sigma_(:,:,tt) + scoreparam(jj)*zeros(k);
        else
            Sigma_(:,:,tt) = Sigma_(:,:,tt) + scoreparam(jj)*ScaledScore(:,:,tt-jj);
        end
    end
    for jj = 1:q
        if (tt-jj) <= 0
            Sigma_(:,:,tt) = Sigma_(:,:,tt) + garchparam(jj)*ini_Sigma;
        else
            Sigma_(:,:,tt) = Sigma_(:,:,tt) + garchparam(jj)*Sigma_(:,:,tt-jj);
        end
    end
    % Random draws
    if strcmp( dist, 'Wish' )
        X(:,:,tt) = matvWishrnd( ...
            Sigma_(:,:,tt)/n, ... % To make Sigma_ the expectation.
            n, ...
            1 ...
        );
        S = score.Sigma_*n;
    elseif strcmp( dist, 'iWish' )
        X(:,:,tt) = matviWishrnd( ...
            Sigma_(:,:,tt)*(n - k - 1), ... % To make Sigma_ the expectation.
            n, ...
            1 ...
        );
        S = score.Sigma_/(n - k - 1);
    elseif strcmp( dist, 'tWish')
        X(:,:,tt) = matvtWishrnd( ...
            Sigma_(:,:,tt)*(nu - 2)/nu/n, ... % To make Sigma_ the expectation.
            n, ...
            nu, ...
            1 ...
        );
        S = score.Sigma_/((nu - 2)/nu/n);
    elseif strcmp( dist, 'itWish')
        X(:,:,tt) = matvitWishrnd( ...
            Sigma_(:,:,tt)*(n-k-1), ... % To make Sigma_ the expectation.
            n, ...
            nu, ...
            1 ...
        );
        S = score.Sigma_/(n-k-1);
    elseif strcmp( dist, 'F' )
        X(:,:,tt) = matvFrnd( ...
            Sigma_(:,:,tt)*(nu - 2)/n, ... % To make Sigma_ the expectation.
            n, ...
            nu, ...
            1 ...
        );
        S = score.Sigma_/((nu - 2)/n);
    elseif strcmp( dist, 'Riesz' )
        cholSig = chol(Sigma_(:,:,tt),'lower');
        Sigma_Riesz = cholSig/diag(n)*cholSig';
        X(:,:,tt) = matvRieszrnd( ...
            Sigma_Riesz, ... % To make Sigma_ the expectation.
            n, ...
            1 ...
        );
        S = score.Sigma*sum(n)/k;
    elseif strcmp( dist, 'iRiesz' )
        invCholInvSig = inv(chol(inv(Sigma_(:,:,tt)),'lower'));
        M = matviRieszexpmat(n);
        Sigma_iRiesz = invCholInvSig'/M*invCholInvSig;
        X(:,:,tt) = matviRieszrnd( ...
            Sigma_iRiesz, ... % To make Sigma_ the expectation.
            n, ...
            1 ...
        );
        S = score.Sigma*sum(diag(M))/k;
    elseif strcmp( dist, 'tRiesz' )
        cholSig = chol(Sigma_(:,:,tt),'lower');
        Sigma_RieszT = cholSig/diag(n)*cholSig'.*(nu - 2)./nu;
        X(:,:,tt) = matvtRieszrnd( ...
            Sigma_RieszT, ... % To make Sigma_ the expectation.
            n, ...
            nu, ...
            1 ...
        );  
        S = score.Sigma*sum(n)/k*nu/(nu - 2);
    elseif strcmp( dist, 'itRiesz' )
        invCholInvSig = inv(chol(inv(Sigma_(:,:,tt)),'lower'));
        M = matviRieszexpmat(n);
        Sigma_iRieszT = invCholInvSig'/M*invCholInvSig;
        X(:,:,tt) = matvitRieszrnd( ...
            Sigma_iRieszT, ... % To make Sigma_ the expectation.
            n, ...
            nu, ...
            1 ...
        );
        S = score.Sigma*sum(diag(M))/k;
    elseif strcmp( dist, 'FRiesz' )
        cholSig = chol(Sigma_(:,:,tt),'lower');
        M = matvFRieszexpmat(n,nu);
        Sigma_FRiesz = cholSig/M*cholSig';
        X(:,:,tt) = matvFRieszrnd( ...
            Sigma_FRiesz, ... % To make Sigma_ the expectation.
            n, ...
            nu, ...
            1 ...
        );
        S = score.Sigma*sum(diag(M))/k;
    end

    % Scaling as in Opschoor et. al (2018)
    S = ivech(S);
        % Opschoor et al. take derivative ignogring that Sigma_ is
        % symmetric, whereas my scores always take symmetry into
        % account. To undo that, half the off-diagonal elements of S.
    S = .5*(ones(p) + eye(p)).*S; 
    ScaledScore(:,:,tt) =  2/(sum(n.*ones(k,1))/k + 1) * Sigma_(:,:,tt)*S*Sigma_(:,:,tt);   

%         % Scaling as proposed by Creal, Koopman, and Lucas (2013)          
%         Score(:,:,tt) = ivech( inv(fisherinfo.Sigma_)*score.Sigma_' );  
        
end
    
end
