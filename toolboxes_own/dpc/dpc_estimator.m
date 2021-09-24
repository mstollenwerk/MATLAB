function [ eparam, logQL, eS, eL, e_d, diagLRL, fcstS, fcstL, fcstd ] = dpc_estimator( covMdata, loadSpec, lP, lQ, compSpec, P, Q, x0 )
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
%   lP           - [OPTIONAL] Non-negative, scalar integer representing the number of covMdata lags in BEKK loading recursion (Default: 1)
%   lQ           - [OPTIONAL] Non-negative, scalar integer representing the number of conditional covariance matix lags in BEKK loading recursion (Default:1)
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
%      [1] Gribisch, B. and Stollenwerk, M. (2016) - Dynamic Principal 
%      Component CAW Models for High-Dimensional Realized Covariance 
%      Matrices. http://ssrn.com/abstract=2888520
%
%      [2] Aielli, G. P. and M. Caporin (2015) - Dynamic Principal
%      Components: a New Class of Multivariate GARCH Models

% Michael Stollenwerk
% michael.stollenwerk@live.com
% 29.12.2016

[n,~,T] = size(covMdata);
%% Defaults
if isempty(lP)
    lP=1;
end
if isempty(lQ)
    lQ=1;
end
if isempty(P)
    P=1;
end
if isempty(Q)
    Q=1;
end
%% x0 - Starting Values for quasi-logLikelihood Optimization
% Create default starting points according to other inputs.
if strcmpi(loadSpec,'scalar') || isempty(loadSpec)
    x0Load = .1*ones(1,lP+lQ);
elseif strcmpi(loadSpec,'restrfull')
    x0Load = ones(1,2*n)*.1;
elseif strcmpi(loadSpec,'ogarch') || strcmpi(loadSpec,'constant')
    x0Load = [0 0];
elseif strcmpi(loadSpec,'rm')
    x0Load = [];
else error('Please give valid loadSpec input: scalar, restrfull, ogarch, constant or []')
end
if strcmpi(compSpec,'garch') || isempty(compSpec)
    x0Comp = ones(1,n*(P+Q))*.01;
elseif strcmpi(compSpec,'har')
    x0Comp = ones(1,n*4)*.1;
else error('Please give valid compSpec input: garch, har or []')
end
% If x0 is provided, check if it has the correct size acc. to other inputs
if ~isempty(x0)
    x0 = x0(n*(n+1)/2+1:end);
    if numel([x0Load, x0Comp])==numel(x0)
        x0Load = x0(1:numel(x0Load));
        x0Comp = x0(1+numel(x0Load):numel(x0Load)+numel(x0Comp));        
    else error('Number of elements in x0 does not fit.')
    end
end
%% Step 1 - Targeting of intercept matrix "Sc"
eSc = mean(covMdata,3);
%% Step 2 - Loadings Recursion Estimation
[ eparamLoad, eL, fcstL ] = dpc_estimatorLoad(covMdata,loadSpec,lP,lQ,x0Load);
%% Step 3 - Components Recursion Estimation
diagLRL = NaN(n,T); % Get component input data
for t=1:T 
    diagLRL(:,t) = diag(eL(:,:,t)'*covMdata(:,:,t)*eL(:,:,t));
end
[ eparamComp, e_d, fcstd ] = dpc_estimatorComp(covMdata,diagLRL,compSpec,P,Q,x0Comp);
%% Output 
eS = NaN(size(covMdata));
for t=1:T % Create conditional covariance matrix "eS" recursion.
    eS(:,:,t) = eL(:,:,t)*diag(e_d(:,t))*eL(:,:,t)';
end
fcstS=fcstL*diag(fcstd)*fcstL';
[ nLogQL, ~ ] = dpcqlike(e_d,diagLRL);
logQL=-nLogQL;
eparam = [vech(chol(eSc,'lower'))' eparamLoad eparamComp];
end
