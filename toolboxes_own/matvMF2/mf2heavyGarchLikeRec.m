function [nLogL,logLcontr,h,z,Corr,H] = mf2heavyGarchLikeRec(param,r,RC)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

m = 63;

[p,~,T] = size(RC);

%% Realized Correlations
sel = tril(true(p),-1);
nCorr = sum(sum(sel));
Rcorr = NaN(nCorr,T);
for tt = 1:T
    rsdevs = sqrt(diag(RC(:,:,tt)));
    CorrMat = RC(:,:,tt)./(rsdevs*rsdevs');
    Rcorr(:,tt) = CorrMat(sel);
end

%% Parameter Transformation
intrcptTau = param(1:p);
aTau = param(p+1);
bTau = param(p+2);
aS = param(p+3);
bS = param(p+4);
aCorr = param(p+5);
bCorr = param(p+6);
% intrcptCorrs = param(p+4+1:p+4+nCorr);
% aCorr = param(p+4+nCorr+1);
% bCorr = param(p+4+nCorr+2);

%% Likelihood Recursion
h = NaN(p,T);
z = NaN(p,T);
v = NaN(p,T);
H = NaN(p,p,T);
Corr = NaN(nCorr,T);

logLcontr = NaN(T,1);

% Initialization
m_ini = ceil(sqrt(T));
w = .06 * .94.^(0:(m_ini-1));
w = w/sum(w);
bachCastRcorr = sum(bsxfun(@times, w, Rcorr(:,1:m_ini)),2);
w = reshape(w,[1 1 m_ini]);
backCastRC = sum(bsxfun(@times,w, RC(:,:,1:m_ini)),3);

h(:,1) = ones(p,1);
z(:,1:m) = repmat(diag(backCastRC),1,m);
v(:,1) = diag(RC(:,:,1))./h(:,1);
H(:,:,1) = backCastRC;
Corr(:,1) = bachCastRcorr;
[~,logLcontr(1)] = mvnormlike(r(:,1), zeros(p,1), H(:,:,1));

for tt = 2:m
    
    h(:,tt) = (1-aS-bS)*ones(p,1) + aS*(diag(RC(:,:,tt-1))./z(:,tt-1)) + bS*h(:,tt-1);
    Dh = diag(sqrt(h(:,tt)));
    
    Dz = diag(sqrt(z(:,tt)));
    
    v(:,tt) = diag(RC(:,:,tt))./h(:,tt);
    
    Corr(:,tt) = (1-aCorr-bCorr)*mean(Rcorr,2) + aCorr*Rcorr(:,tt) + bCorr*Corr(:,tt-1);
    
    H_ = .5.*eye(p);
    H_(sel) = Corr(:,tt);
    H_ = H_ + H_';
    H(:,:,tt) = Dz*Dh*H_*Dh*Dz;
    
    try
        [~,logLcontr(tt)] = mvnormlike(r(:,tt), zeros(p,1), H(:,:,tt));
    catch
        nLogL = inf;
        logLcontr = NaN;
        h = NaN;
        z = NaN;
        Corr = NaN;
        H = NaN;
        return
    end
    
end

for tt = m+1:T
    
    h(:,tt) = (1-aS-bS)*ones(p,1) + aS*(diag(RC(:,:,tt-1))./z(:,tt-1)) + bS*h(:,tt-1);
    Dh = diag(sqrt(h(:,tt)));
    
    v(:,tt) = diag(RC(:,:,tt))./h(:,tt);
    
    z(:,tt) = intrcptTau + aTau*mean(v(:,tt-m:tt-1),2) + bTau*z(:,tt-1);
    Dz = diag(sqrt(z(:,tt)));
    
    Corr(:,tt) = (1-aCorr-bCorr)*mean(Rcorr,2) + aCorr*Rcorr(:,tt) + bCorr*Corr(:,tt-1);
    
    H_ = .5.*eye(p);
    H_(sel) = Corr(:,tt);
    H_ = H_ + H_';
    H(:,:,tt) = Dz*Dh*H_*Dh*Dz;
    
    try
        [~,logLcontr(tt)] = mvnormlike(r(:,tt), zeros(p,1), H(:,:,tt));
    catch
        nLogL = inf;
        logLcontr = NaN;
        h = NaN;
        z = NaN;
        Corr = NaN;
        H = NaN;
        return
    end
    
end

% Loglikelihood
nLogL = -sum(logLcontr);

end