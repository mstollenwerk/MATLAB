function [nLogL,logLcontr,h,tau] = mf2heavyGarchLikeRecStep1(param,r,RC,m)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

[p,~,T] = size(RC);
%% Parameter Transformation
intrcptTau = param(1:p);
aTau = param(p+1);
bTau = param(p+2);
aS = param(p+3);
bS = param(p+4);

%% Likelihood Recursion
h = NaN(p,T);
tau = NaN(p,T);
v = NaN(p,T);

logLcontr = NaN(T,1);

% Initialization
m_ini = ceil(sqrt(T));
w = .06 * .94.^(0:(m_ini-1));
w = w/sum(w);
w = reshape(w,[1 1 m_ini]);
backCastRC = sum(bsxfun(@times,w, RC(:,:,1:m_ini)),3);

h(:,1) = ones(p,1);
% tau(:,1:m) = repmat(var(r,[],2),1,m);
tau(:,1:m) = repmat(diag(backCastRC),1,m);
v(:,1) = diag(RC(:,:,1))./h(:,1);
sig_sq = h(:,1).*tau(:,1);
logLcontr(1) = -.5*sum(log(2*pi) + log(sig_sq) + (r(:,1).^2)./sig_sq);

for tt = 2:T
    
    h(:,tt) = (1-aS-bS)*ones(p,1) + aS*(diag(RC(:,:,tt-1))./tau(:,tt-1)) + bS*h(:,tt-1);
    v(:,tt) = diag(RC(:,:,tt))./h(:,tt);
    if tt > m
        tau(:,tt) = intrcptTau + aTau*mean(v(:,tt-m:tt-1),2) + bTau*tau(:,tt-1);
    end
    sig_sq = h(:,tt).*tau(:,tt);
    
    try
        logLcontr(tt) = -.5*sum(log(2*pi) + log(sig_sq) + (r(:,tt).^2)./sig_sq);
    catch
        nLogL = inf;
        logLcontr = NaN;
        h = NaN;
        tau = NaN;
        return
    end
    
end

% Loglikelihood
nLogL = -sum(logLcontr);

end

function out_ = decay(omega,L,m,X,varargin)
% The decay function takes as arguments the parameters omega, 
%"numer of lags" L, "averaging length" m, "data to average" X and in
% case L*m is larger than length of data you can put the optional argument
% how you want to initialize X.

length_ = size(X,2);
if nargin == 5
    X = cat(2,repmat(varargin{1},1,L*m-length_),X);
end

if size(X,2) ~= L*m
    error('Something went wrong.')
end
X = flip(X,2);

rho = NaN(L,1);
if L == 1
    rho = 1;
else
    for l = 1:L
        %weights
        rho(l) = (1-l/L)^(omega-1);
    end
    rho = rho./sum(rho);
end
%rvm
out_ = 0;
for l = 1:L
    out_ = out_ + rho(l)*mean(X(:,m*(l-1)+1:m*l),2);
end

end
