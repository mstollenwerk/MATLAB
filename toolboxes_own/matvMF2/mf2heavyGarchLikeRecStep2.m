function [nLogL,logLcontr,RCorrE] = mf2heavyGarchLikeRecStep2(param,rStan,RCorr,m)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

omega = 1;
L = 1;

[p,~,T] = size(RCorr);

%% Parameter Transformation
aCorr = param(1);
bCorr = param(2);

%% Likelihood Recursion
RCorrE = NaN(p,p,T);
logLcontr = NaN(T,1);

% Initialization
m_ini = ceil(sqrt(T));
w = .06 * .94.^(0:(m_ini-1));
w = w/sum(w);
w = reshape(w,[1 1 m_ini]);
backCastRCorr = sum(bsxfun(@times,w, RCorr(:,:,1:m_ini)),3);

RCorrE(:,:,1) = backCastRCorr;
logLcontr(1) = ...
            - logdet(RCorrE(:,:,1)) ...
            - rStan(:,1)'*inv(RCorrE(:,:,1))*rStan(:,1) ...
            + rStan(:,1)'*rStan(:,1);
        
for tt = 2:T
   
    RCorrE(:,:,tt) = (1-aCorr-bCorr)*decay(omega,L,m,RCorr(:,:,max(tt-L*m,1):tt-1),backCastRCorr) ...
                        + aCorr*RCorr(:,:,tt-1) + bCorr*RCorrE(:,:,tt-1);
    
    try
        logLcontr(tt) = ...
            - logdet(RCorrE(:,:,tt)) ...
            - rStan(:,tt)'*inv(RCorrE(:,:,tt))*rStan(:,tt) ...
            + rStan(:,tt)'*rStan(:,tt);
    catch
        nLogL = inf;
        logLcontr = NaN;
        RCorrE = NaN;
        return
    end
    
end

% Loglikelihood
nLogL = -.5*sum(logLcontr);

end

function out_ = decay(omega,L,m,X,varargin)
% The decay function takes as arguments the parameters omega, 
%"numer of lags" L, "averaging length" m, "data to average" X and in
% case L*m is larger than length of data you can put the optional argument
% how you want to initialize X.

length_ = size(X,3);
if nargin == 5 && length_ < L*m
    X = cat(3,repmat(varargin{1},1,1,L*m-length_),X);
end

if size(X,3) ~= L*m
    error('Something went wrong.')
end
X = flip(X,3);

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
    out_ = out_ + rho(l)*mean(X(:,:,m*(l-1)+1:m*l),3);
end

end