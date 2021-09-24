function [ eparamComp, eS, eL, e_d, diagLRL, fcstS, fcstL, fcstd ] = dpc_estimator_Q( covMdata, eQ, fcstQ, P, Q, factors )
%DPCCAWESTIMATOR estimates the DPC-CAW model using a 3-step estimation.
%   DPCCAWESTIMATOR 
%
% USAGE:
%  [EPARAM] = dpccawestimator(COVMDATA,[],P,O,Q)
%  [PARAMETERS,LogQL,eS,eL,e_d] = dpccawestimator(DATA,DATAASYM,P,O,Q,TYPE,STARTINGVALS,OPTIONS)
%
% INPUTS:
%   COVMDATA     - N by N by T array of covariance estimators (e.g. realized covariance)
%   LOADSPEC     - [OPTIONAL] String, one of 'Scalar' (Default) ,'RestrFull', 'ogarch' or 'constant'
%   COMPSPEC     - [OPTIONAL] String, one of 'garch' (Default) or 'har'
%   P            - [OPTIONAL] Non-negative, scalar integer representing the number of innovations in components recursion (Default: 1)
%   Q            - [OPTIONAL] Non-negative, scalar integer representing the number of conditional variance lags in components recursion (Default:1)
%   X0           - [OPTIONAL] Starting values for optimization. Format see
%                  EPARAM in OUTPUTS
%
% OUTPUTS:
%   EPARAM       - Vector of parameters. The form of the parameters depends on LOADSPEC and COMPSPEC.
%                  Example:                    
%                    'Scalar' and 'garch':
%                    [vech(LoadIntercept)' a_load b_load a(1) ... a(p) b(1) ... b(q)]'  (all scalars)
%   LogQL        - The quasi log likelihood at the optimum
%   eS           - [N N T] dimension matrix of conditional covariances
%   eL           - [N N T] dimension matrix of conditional loadings
%   e_d          - [N T] dimension matrix of conditional components
%
% COMMENTS:
%   The dynamics are given by
%
%  See also 
%
% REFERENCES:
%      [1] Aielli, G. P. and M. Caporin (2015) - Dynamic Principal
%      Components: a New Class of Multivariate GARCH Models

% Michael Stollenwerk
% michael.stollenwerk@live.com
% 28.08.2016

[n,~,T] = size(covMdata);
%% Defaults
if isempty(P)
    P=1;
end
if isempty(Q)
    Q=1;
end

[eL_,~]=dpceig(mean(covMdata,3)); % Get Intercept Eigenvectors to enforce sign convention below. (Assumptions 2.1.2. - 2.1.3.)
eL = NaN(size(eQ)); % Storage
for t=1:size(eQ,3)
    [eL(:,:,t),~] = dpceig(eQ(:,:,t),eL_); 
end
[fcstL,~] = dpceig(fcstQ,eL_);

%% Components Recursion Estimation
diagLRL = NaN(n,T); % Get component input data
for t=1:T 
    diagLRL(:,t) = diag(eL(:,:,t)'*covMdata(:,:,t)*eL(:,:,t));
end

if isempty(factors)
    [ eparamComp, e_d, fcstd ] = dpc_estimatorComp(covMdata,diagLRL,[],P,Q,[]);
elseif factors==0
    [ eparamComp, e_d, fcstd ] = dpc_estimatorComp_restr( covMdata, diagLRL, [], P, Q, [], factors );
else error('Only [] and 0 accepted as input for factors')
end
%% Output 
eS = NaN(size(covMdata));
for t=1:T % Create conditional covariance matrix "eS" recursion.
    eS(:,:,t) = eL(:,:,t)*diag(e_d(:,t))*eL(:,:,t)';
end
fcstS=fcstL*diag(fcstd)*fcstL';
end
