function [ nLogL, h, tau, V_m ] = mf2_garch_likeRec(parameters, rv)
    
%     M = 5;
%     L = 12;
%     m = L*M; 
    m = 63;

    alpha = parameters(1);
    beta = parameters(2);
    
    lambda_0 = parameters(3);
    lambda_1 = parameters(4);
    lambda_2 = parameters(5);
    
    n = parameters(6);

    h = ones(size(rv));
    
    tau = ones(size(rv))*mean(rv);
    V = zeros(size(rv));
    V_m = zeros(size(rv));
    
    for t = 2:m
        
       h(t) = (1-alpha-beta) + alpha .* rv(t-1) ./ tau(t-1) + beta .* h(t-1); 
       
    end
    
    
    for t = (m+1):length(rv)

       h(t) = (1-alpha-beta) + alpha .* rv(t-1) ./ tau(t-1) + beta .* h(t-1); 
       
       V(t) = rv(t) ./ h(t);
       V_m(t) = sum(V(t-(m-1):t))./m;
       tau(t) = lambda_0 + lambda_1 * V_m(t-1) + lambda_2 * tau(t-1);
%         tau(t) = lambda_0 + lambda_1 * decay(lambda_2,L,M,flip(rv(1:t-1))); 
       
    end
    
    for t = 1:length(rv)
        try
            lls(t) = matvsWishlike(h(t)*tau(t),n,rv(t));
        catch
            nLogL = inf;
            h = NaN;
            tau = NaN;
            V_m =NaN;
            return
        end
    end
    nLogL = sum(lls);
    
end


function out_ = decay(omega,L,m,rv)

rho = NaN(L,1);
for l = 1:L
    %weights
    rho(l) = (1-l/L)^(omega-1);
end
rho = rho./sum(rho);
%rvm
out_ = 0;
for l = 1:L
    out_ = out_ + rho(l)*mean(rv(m*(l-1)+1:m*l));
end

end

function [ nLogL, logLcontr ] = matvsWishlike( Sigma_, n, X )
%MATVSWISHLIKE Log-likelihood of Wishart distribution.

[p,~,N] = size(X);
%% Log-Likelihood
logLcontr = NaN(N,1);

log_norm_const = ...
    n*p/2*log(n) ...
  - n*p/2*log(2) ...
  - mvgammaln(n/2, p) ...
  - n/2*log(det(Sigma_));
	  
for ii = 1:N
    R = X(:,:,ii);
 
    log_kernel = ...
      -  trace(n*(Sigma_\R)) ...
      + (n-p-1)*log(det(R));
		
    logLcontr(ii) = log_norm_const + .5*log_kernel;
end
nLogL = -sum(logLcontr);
end

function y = mvgammaln(x,d)
%MVGAMMALC computes the natural logarithm of the multivariate gamma function 
%
% USAGE:
%  Y = mvgammaln(x,d)
%
% INPUTS:
%   X            - x-value
%   D            - "number of sums" parameter
%
% OUTPUTS:
%   Y            - y-value
%
% COMMENTS:
%   Used in the probability density function of the Wishart and inverse 
%   Wishart distributions.
%   Gamma_d(x) = pi^(d(d-1)/4) \prod_(j=1)^d Gamma(x+(1-j)/2)
%   log(Gamma_d(x)) = d(d-1)/4 log(pi) + \sum_(j=1)^d log(Gamma(x+(1-j)/2))
%
% REFERENCES:
%      [1] James (1964) - Distributions of Matrix Variates and Latent 
%      Roots Derived from Normal Samples. 

% Michael Stollenwerk
% michael.stollenwerk@live.com
% 12.02.2017

y = d*(d-1)/4*log(pi)+sum(gammaln(x+(1-(1:d))/2));

end