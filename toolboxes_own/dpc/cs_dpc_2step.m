function [ eM, eparam, nLogL, logLcontr, e_expData, eS, eQ, eL, e_d, e_diagLRL, fcstData ] =...
    cs_dpc_2step( data, Pl, Ql, Pc, Qc, x0 )
%DPCCAW estimates the consistent DPC-CAW model.
%
% USAGE:
%       [ eM, eparam, logL, eS, eQ, e_d, fcstS, enu ] = cs_dpc( data, Pl, Ql, Pc, Qc, x0 )
%
% INPUTS:
%       data       - Data Input, i.e. realized cov-matrices
%       Pl         - Number of lagged covMats in Loadings Recursion
%       Ql         - Number of lagged Q Matrices in Loadings Recursion
%       Pc         - Number of lagged diagLRL in Component Recursion
%       Qc         - Number of lagged d_t in Component Recursion
%       x0         - [Optional] Starting point for the minimization
%
% OUTPUTS:
%       eM         - estimated M matrix of consistent DPC-CAW (see paper)
%       eparam     - estimated parameters of centered process
%       nLogL      - value of negative quasi log-likelihood
%       logLcontr  - log-likelihood contributions
%       e_expData  - expected value of data, scale matrix of wishart
%       eS         - estimated expectations for Cov-Matrices "S"
%       eQ         - estimated loadings recursion matrices
%       eL         - eigenvalue matrix recursion
%       e_d        - estimated components recursion
%       e_diagLRL  - estimated components data recursion
%       fcstData   - Forecast of Data in T+1
%       enu        - estimated dof parameter of centered process
%
% COMMENTS:
%
%  See also 
%
% REFERENCES:
%   [1] Gribisch and Stollenwerk

% Michael Stollenwerk
% michael.stollenwerk@live.com
% 12.06.2017
% 08.08.2017: removed helper functions
%             included eM into likelihood

[n,~,~] = size(data);
%% Step 1 - center data
eM = mean(data,3);
%% x0
if isempty(x0)
    x0 = [ones(1,Pl)*.05/Pl ones(1,Ql)*.92/Ql ones(1,Pc*n)*.05/Pc ones(1,Qc*n)*.92/Qc];
end
%% Minimization 
% Restrictions for Loadings Parameters
if Pl>1 || Ql>1
    error('Only P=Q=1 supported in the loadings recursion so far, because I dunno the restrictions for bigger P/Q.')
end

AA = [1 1 zeros(1,n*(Pc+Qc));
      zeros(n,2) eye(n) eye(n)]; % Stationarity of Loadings Recursion (a+b <=1) 
bb = ones(1+n,1); %*(1-1e-8);
lb = zeros(Pl+Ql+n*(Pc+Qc),1);
ub = ones(Pl+Ql+n*(Pc+Qc),1);

options = optimoptions('patternsearch','Display','iter'); %,'StepTolerance',1e-5 ,'MaxFunctionEvaluations',2e3
eparam = patternsearch(@(param)cs_dpcLikeRec(param,eM,data,Pl,Ql,Pc,Qc), x0,...
    AA, bb,[],[],lb,ub,[], options);

% Optimization w.r.t. dof "nu".
% options = optimoptions('fmincon','Display','iter');
[nLogL, logLcontr, e_expData, eS, eQ, eL, e_d, e_diagLRL, fcstData ] =...
    cs_dpcLikeRec( eparam, eM, data, Pl, Ql, Pc, Qc);
% [enu, nLogL] = fmincon(@(nu)cawlike_givn2ndpart_( nu, data_, LogLcontr_2ndpart ), 2*n,[],[],[],[],n-1,[],[], options);
% logL = -nLogL;
end
