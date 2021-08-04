function [ nLogL, logLcontr ] = mLaplacelike( data_mat, Sigma_, Phi_ )
%MLAPLACELIKE
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
%      [1]                     
%
% Michael Stollenwerk
% michael.stollenwerk@live.com
% 20.02.2020
%
% DEPENDENCIES:
%
%%
[p,n,K] = size(data_mat);
[p_,~,K_] = size(Sigma_);
[n_,~,K__] = size(Phi_);

if p~=p_ || K~=K_ || K_~=K__ || n~=n_
    error('Input array dimensions do not agree.')
end

logLcontr = NaN(K,1);
for kk = 1:K
    X = data_mat(:,:,kk);
    c = log(2) - p*n/2*log(2*pi) - n/2*log(det(Sigma_(:,:,kk))) - ...
        p/2*log(det(Phi_(:,:,kk)));
    trQ = trace((Sigma_(:,:,kk)\X)*(Phi_(:,:,kk)\X'));
    log_kernel = (2-p*n)/4*log(trQ/2) + log(besselk((2-p*n)/2,sqrt(2*trQ)));
    
    logLcontr(kk) = c + log_kernel;
end
nLogL = -sum(logLcontr);

% vecX = NaN(K,p*n);
% Sigma_vecX = NaN(p*n,p*n,K);
% for ii = 1:K
%     X = data_mat(:,:,ii);
%     vecX(ii,:) = X(:);
%     Sigma_vecX(:,:,ii) = kron(Sigma_(:,:,ii), Phi_(:,:,ii));
% end
% [nLogL, logLcontr] = vLaplacelike(vecX, Sigma_vecX);
end
