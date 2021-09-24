function [ S, L, d, diagLRL, fcstS, fcstL, fcstd ] = dpc_estimator_fcst( covMdata, paramLoad, paramComp, P, Q, num_comp, type )
%DPCCAWESTIMATOR estimates the DPC-CAW model using a 3-step estimation.
%   DPCCAWESTIMATOR 
%
% USAGE:
%  [S,L,d,diagLRL,fcstS,fcstL,fcstd] = dpc_estimator_fcst(COVMDATA,PARAMLOAD,PARAMCOMP,P,Q,NUM_COMP,TYPE)
%
% INPUTS:
%   COVMDATA     - N by N by T array of covariance estimators (e.g. realized covariance)
%   PARAMLOAD    - Loading Recursion parameters
%   PARAMCOMP    - Component Recursions parameters
%   P            - Non-negative, scalar integer representing the number of innovations in components recursion
%   Q            - Non-negative, scalar integer representing the number of conditional variance lags in components recursion
%   NUM_COMP     - Number of components being modeled with GARCH specification
%   TYPE         - Describes how the non-GARCH modeled components are determined either they are set to unconditional means or S and fcstS
%                  are created from NUM_COMP components only
%
% OUTPUTS:
%
% COMMENTS:
%   The dynamics are given by
%
%  See also 
%
% REFERENCES:
%      [1] Gribisch, B. and Stollenwerk, M. (2016) - Dynamic Principal 
%      Component CAW Models for High-Dimensional Realized Covariance 
%      Matrices
%      [2] Aielli, G. P. and M. Caporin (2015) - Dynamic Principal
%      Components: a New Class of Multivariate GARCH Models

% Michael Stollenwerk
% michael.stollenwerk@live.com
% 17.12.2016

[N,~,T] = size(covMdata);
if size(paramLoad,1)>size(paramLoad,2)
    paramLoad=paramLoad';
end
if size(paramComp,1)>size(paramComp,2)
    paramComp=paramComp';
end
%% Step 1 - Targeting of intercept matrix "Sc"
eSc = mean(covMdata,3);
%% Step 2 - Get Loadings Recursion and Forecast
[eQ, fcstQ] = sbekkRec([vech(chol(eSc.*(1-sum(paramLoad)),'lower'))', paramLoad],covMdata,1,1);
eQ=reshape(cell2mat(eQ),N,N,T);
%% Output
[eL_,edc]=dpceig(mean(covMdata,3)); % Get Intercept Eigenvectors to enforce sign convention below. (Assumptions 2.1.2. - 2.1.3.)
L = NaN(size(eQ)); % Storage
for t=1:size(eQ,3)
    [L(:,:,t),~] = dpceig(eQ(:,:,t),eL_); 
end
[fcstL,~] = dpceig(fcstQ,eL_);
%% Step 3 - Get Component Recursions and Forecasts
diagLRL = NaN(N,T); % Get component input data
for t=1:T 
    diagLRL(:,t) = diag(L(:,:,t)'*covMdata(:,:,t)*L(:,:,t));
end
%% Storage and Reshape x0 for individual deployment
d = repmat(edc,1,T);
fcstd=edc;
paramComp = reshape(paramComp,N,[]); % Each paramComp row now corresponds to component parameters for one component.
for i=1:num_comp
    [d(i,:), fcstd(i)] = garchRec([(1-sum(paramComp(i,:)))*edc(i) paramComp(i,:)],diagLRL(i,:),P,Q);
end
%% Output 
S = NaN(size(covMdata));
if strcmpi(type,'const')
    for t=1:T % Create conditional covariance matrix "eS" recursion.
        S(:,:,t) = L(:,:,t)*diag(d(:,t))*L(:,:,t)';
    end
    fcstS=fcstL*diag(fcstd)*fcstL';
elseif strcmpi(type,'none')
    for t=1:T % Create conditional covariance matrix "eS" recursion.
        S(:,:,t) = L(:,1:num_comp,t)*diag(d(1:num_comp,t))*L(:,1:num_comp,t)';
    end
    fcstS=fcstL(:,1:num_comp)*diag(fcstd(1:num_comp))*fcstL(:,1:num_comp)';
else error('Specify valid type of restriction: const or none')
end

end
