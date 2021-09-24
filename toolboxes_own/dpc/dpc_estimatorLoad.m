function [ eparamLoad, eL, fcstL] = dpc_estimatorLoad( covMdata, loadSpec , lP, lQ, x0Load )
%DPCCAWESTIMATORLOAD estimates loadings recursion in DPCCAW
%   DPCCAWESTIMATORLOAD
%
% USAGE:
%  [EPARAMLOAD] = 
%  [EPARAMLOAD,EL,FCSTL] = 
%
% INPUTS:
%   COVMDATA     - N by N by T array of covariance estimators (e.g. realized covariance)
%   LOADSPEC     - [OPTIONAL] String, one of 'Scalar' (Default) ,'RestrFull', 'ogarch' or 'constant'
%
% OUTPUTS:
%   eL           - [N N T] dimension matrix of conditional loadings
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

%% Optimizations
if strcmpi(loadSpec,'scalar') || isempty(loadSpec)
    [ eparamLoad, ~, eQ, fcstQ ] = sbekk_targ( covMdata, [], lP, lQ, x0Load ); % Estimate stationary scalar CAW (a.k.a. stationary scalar Re-BEKK) with covariance targeting.
% elseif strcmpi(loadSpec,'ogarch') || strcmpi(loadSpec,'constant')
%     eparamLoad = [0 0];
%     eQ = repmat(mean(covMdata,3),1,1,size(covMdata,3)); %such that in the main function 'dpccawestimator' eL is created correctly.
elseif strcmpi(loadSpec,'restrfull') % COMPOSITE LIKELIHOOD METHOD USED
    [n,~,T] = size(covMdata);
%     AA = [-eye(2*n) % all paramters >= 0 except Sc; actually at least one a_i and b_i is required to be strictly >0, but ok... see HVMA p.51
%           eye(n) eye(n)]; % stationarity (a_i+b_i <=0) Assumption 2.6.1 DPC paper Aielli and Caporin (2016)
%     bb = [zeros(1,2*n) ones(1,n)];
    RCdata = mat2cell(covMdata(:,:),n,ones(1,T)*n);
    Sc = mean(covMdata,3);
    options = optimoptions('fminunc','Display','iter-detailed','UseParallel','always'); % ,'MaxFunctionEvaluations'
    eparamLoad = fminunc(@(paramLoad)cawcompqlike_full(bekkRecVechC11(parbekk2rf([vech(Sc)' paramLoad],n),RCdata),RCdata,'diagonal'), x0Load, options);     
    [eQ, fcstQ] = bekkRecVechC11(parbekk2rf([vech(Sc)' eparamLoad],n),RCdata);
    eQ = reshape(cell2mat(eQ),n,n,T);
elseif strcmpi(loadSpec,'rm')
    tau0 = 1560;         % Does not matter
    rho  = 1;            % Does not matter
    tau1 = -1/log(.94);
    kmax = 1;
    eQ = riskmetrics2006(covMdata,tau0,tau1,kmax,rho);
    fcstQ = .94*eQ(:,:,end) +.06*covMdata(:,:,end);
    eparamLoad=[];
else error('Invalid loading specification')
end
%% Output
[eL_,~]=dpceig(mean(covMdata,3)); % Get Intercept Eigenvectors to enforce sign convention below. (Assumptions 2.1.2. - 2.1.3.)
eL = NaN(size(eQ)); % Storage
for t=1:size(eQ,3)
    [eL(:,:,t),~] = dpceig(eQ(:,:,t),eL_); 
end
[fcstL,~] = dpceig(fcstQ,eL_);
end
    